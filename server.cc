/*
    g++ server.cc -o server -lpthread -I ./tfhe/src/include/ -L ./tfhe/build/libtfhe/ -ltfhe-spqlios-fma -O3
 */

#include <cstdio>
#include <cassert>
#include <unistd.h>
#include <vector>
#include <queue>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <memory>
#include <atomic>
#include <tfhe.h>
#include <tfhe_io.h>

#include "constants.h"
using namespace std;

using Size = uint32_t;

class BlockingQueue {
public:
  BlockingQueue() : shutdown_(false) {}

  void Push(Size v) {
    bool was_empty = false;

    {
      lock_guard<mutex> lock(mtx_);
      was_empty = queue_.empty();
      queue_.push(v);
    }

    if (was_empty) {
      cv_.notify_one();
    }
  }

  bool Pop(Size &ref) {
    std::unique_lock<std::mutex> lock(mtx_);

    while (!shutdown_ && queue_.empty()) {
      cv_.wait(lock);
    }

    if (shutdown_) return false;

    ref = queue_.front();
    queue_.pop();
    return true;
  }

  void Shutdown() {
    shutdown_ = true;
    cv_.notify_all();
  }

private:
  std::mutex mtx_;
  std::condition_variable cv_;
  std::queue<Size> queue_;
  std::atomic<bool> shutdown_;
};

const char g_num_fanin[18] = {-1, 0, 0, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3};

class Executor {
private: 
  Size N_;
  Size Q_;
  Size O_;
  vector<Size> indices_out_;
  vector<LogicType> logic_type_;
  vector<vector<Size>> es_;
  vector<array<Size, 3>> res_;
  vector<shared_ptr<atomic<char>>> num_precond_;
  LweSample *wires_;

  TFheGateBootstrappingCloudKeySet* bk_;

  void ParseDatabase(char *database_path) {
    FILE *fp = fopen(database_path, "rb");
    if (!fp) {
      printf("Could not open the database file %s\n", database_path);
      exit(-1);
    }

    auto u8 = [&fp, database_path]() {
      uint8_t val;
      int res = fread(&val, sizeof(uint8_t), 1, fp);
      if (res != 1) {
        printf("The database file %s is corrupted\n", database_path);
        exit(-1);
      }

      return val;
    };

    auto u16 = [&fp, database_path]() {
      uint16_t val;
      int res = fread(&val, sizeof(uint16_t), 1, fp);
      if (res != 1) {
        printf("The database file %s is corrupted\n", database_path);
        exit(-1);
      }

      return val;
    };

    auto u32 = [&fp, database_path]() {
      uint32_t val;
      int res = fread(&val, sizeof(uint32_t), 1, fp);
      if (res != 1) {
        printf("The database file %s is corrupted\n", database_path);
        exit(-1);
      }

      return val;
    };

    auto u64 = [&fp, database_path]() {
      uint64_t val;
      int res = fread(&val, sizeof(uint64_t), 1, fp);
      if (res != 1) {
        printf("The database file %s is corrupted\n", database_path);
        exit(-1);
      }

      return val;
    };

    auto magic = u32();
    if (magic != 0xdeadbeefUL) {
      printf("The file %s is not a database file or was saved in the other endian\n", database_path);
      exit(-1);
    }

    Q_ = u32();
    if (Q_ < 0) {
      printf("Give me a break. Q_ < 0...\nThe file %s is broken, or you gave too long strings...\n", database_path);
      exit(-1);
    }

    assert(g_num_fanin[INPUT] == 0);
    assert(g_num_fanin[CONSTANT] == 0);
    assert(g_num_fanin[NOT] == 1);
    assert(g_num_fanin[AND] == 2);
    assert(g_num_fanin[OR] == 2);
    assert(g_num_fanin[XOR] == 2);
    assert(g_num_fanin[NAND] == 2);
    assert(g_num_fanin[NOR] == 2);
    assert(g_num_fanin[XNOR] == 2);
    assert(g_num_fanin[ANDNY] == 2);
    assert(g_num_fanin[ANDYN] == 2);
    assert(g_num_fanin[ORNY] == 2);
    assert(g_num_fanin[ORYN] == 2);
    assert(g_num_fanin[MUX] == 3);

    N_ = 0;
    Size nmax = 0;
    logic_type_.reserve(Q_);
    num_precond_.reserve(Q_);
    es_.reserve(Q_);
    res_.reserve(Q_);
    while (1) {
      LogicType lt = (LogicType)u8();
      if (lt == ERROR) break;

      logic_type_.emplace_back(lt);
      char sz = g_num_fanin[logic_type_[N_]];
      num_precond_.emplace_back(shared_ptr<atomic<char>>(new atomic<char>(sz)));

      array<Size, 3> res;
      for (auto j=0; j<sz; j++) {
        Size v = u32();
        res[j] = v;

        if (nmax < v+1) {
          nmax = v+1;
          es_.resize(nmax);
        }
        es_[v].emplace_back(N_);
      }
      res_.emplace_back(move(res));
      ++N_;
    }
    fclose(fp);
  }

  void ParseIndices(char *indices_path) {
    FILE *fp = fopen(indices_path, "rb");
    if (!fp) {
      printf("Could not open the indices file %s\n", indices_path);
      exit(-1);
    }

    auto u32 = [&fp, indices_path]() {
      uint32_t val;
      int res = fread(&val, sizeof(uint32_t), 1, fp);
      if (res != 1) {
        printf("The indices file %s is corrupted\n", indices_path);
        exit(-1);
      }

      return val;
    };

    indices_out_.reserve(Q_);
    O_ = 0;
    while (1) {
      auto idx = u32();
      if (idx == 0) break;
      indices_out_.emplace_back(idx-1);
      ++O_;
    }
  }

public:
  Executor(char *key_path, char *query_path, char *database_path, char *indices_path) {
    ParseDatabase(database_path);
    ParseIndices(indices_path);

    FILE* cloud_key = fopen(key_path, "rb");
    if (!cloud_key) {
      printf("Could not open the cloud key file %s\n", key_path);
      exit(-1);
    }

    bk_ = new_tfheGateBootstrappingCloudKeySet_fromFile(cloud_key);
    if (!bk_) {
      printf("Could not prepare the cloud key from the file %s\n", key_path);
      exit(-1);
    }
    fclose(cloud_key);

    FILE *query_fp = fopen(query_path, "rb");
    if (!query_fp) {
      printf("Could not open the encrypted query file %s\n", query_path);
      exit(-1);
    }
    
    wires_ = new_gate_bootstrapping_ciphertext_array(N_, bk_->params);
    if (!wires_) {
      printf("Could not allocate memory for wires_\n");
      exit(-1);
    }

    for (Size i=0; i<Q_; i++) {
      import_gate_bootstrapping_ciphertext_fromFile(query_fp, &wires_[i], bk_->params);
    }
    fclose(query_fp);
  }

  void Eval(char *result_path) {
    atomic<Size> num_remain(N_);
    BlockingQueue que;

    for (Size i=0; i<Q_; i++) {
      assert(logic_type_[i] == INPUT);
      --num_remain;
      for (Size v : es_[i]) {
        if (--(*num_precond_[v]) == 0) {
          que.Push(v);
        }
      }
    }

    for (Size i=Q_; i<Q_+2; i++) {
      assert(logic_type_[i] == CONSTANT);
      bootsCONSTANT(&wires_[i], i-Q_, bk_);
      --num_remain;
      for (Size v : es_[i]) {
        if (--(*num_precond_[v]) == 0) {
          que.Push(v);
        }
      }
    }

    auto num_thread = thread::hardware_concurrency();
    num_thread = num_thread * 2.0;
    vector<thread> threads(num_thread);
    for (auto i=0; i<num_thread; i++) {
      threads[i] = thread(
        [this, &num_remain, &que] {
          while (1) {
            Size v;
            bool got_val = que.Pop(v);
            if (!got_val) return;
            
            --num_remain;

            auto &args = res_[v];

            switch (logic_type_[v]) {
            case NOT:
              bootsNOT(&wires_[v], &wires_[args[0]], bk_);
              break;

            case AND:
              bootsAND(&wires_[v], &wires_[args[0]], &wires_[args[1]], bk_);
              break;

            case OR:
              bootsOR(&wires_[v], &wires_[args[0]], &wires_[args[1]], bk_);
              break;

            case XOR:
              bootsXOR(&wires_[v], &wires_[args[0]], &wires_[args[1]], bk_);
              break;

            case NAND:
              bootsNAND(&wires_[v], &wires_[args[0]], &wires_[args[1]], bk_);
              break;

            case NOR:
              bootsNOR(&wires_[v], &wires_[args[0]], &wires_[args[1]], bk_);
              break;

            case XNOR:
              bootsXNOR(&wires_[v], &wires_[args[0]], &wires_[args[1]], bk_);
              break;

            case ANDNY:
              bootsANDNY(&wires_[v], &wires_[args[0]], &wires_[args[1]], bk_);
              break;

            case ANDYN:
              bootsANDYN(&wires_[v], &wires_[args[0]], &wires_[args[1]], bk_);
              break;

            case ORNY:
              bootsORNY(&wires_[v], &wires_[args[0]], &wires_[args[1]], bk_);
              break;

            case ORYN:
              bootsORYN(&wires_[v], &wires_[args[0]], &wires_[args[1]], bk_);
              break;

            case MUX:
              bootsMUX(&wires_[v], &wires_[args[0]], &wires_[args[1]], &wires_[args[2]], bk_);
              break;

            default:
              assert(0);
            }

            if (es_.size() <= v) continue;
            for (Size u : es_[v]) {
              if (--(*num_precond_[u]) == 0) {
                que.Push(u);
              }
            }
          }
        }
      );
    }

    while (num_remain > 0) {
      printf("%u evaluations remaining...\n",  num_remain.load());
      sleep(10);
    }

    printf("Done! Waiting for threads.join...\n");
    que.Shutdown();
    for (auto i=0; i<num_thread; i++) {
      threads[i].join();
    }

    printf("Done! Writing out the results...\n");

    FILE *result_fp = fopen(result_path, "wb");
    if (!result_fp) {
      printf("Could not save the answer to %s.\n", result_path);
      exit(-1);
    }

    for (auto i=0; i<O_; i++) {
      if (i%100 == 0) printf("%u bits remaining ...\n", O_-i);
      Size idx = indices_out_[i];
      export_gate_bootstrapping_ciphertext_toFile(result_fp, &wires_[idx], bk_->params);
    }
    fclose(result_fp);
  }

  ~Executor() {
    delete_gate_bootstrapping_ciphertext_array(N_, wires_);
    delete_gate_bootstrapping_cloud_keyset(bk_);
  }
};

int main(int argc, char *argv[]) {
  setbuf(stdout, NULL);

  if (argc < 6) {
    assert(argc >= 1);
    printf("Usage: %s (cloud_key_path) (enc_query_path) (database_path) (indices_path) (path_for_enc_answer)\n", argv[0]);
    exit(0);
  }

  printf("Constructing a logic circuit from the database file...\n");
  Executor executor(argv[1], argv[2], argv[3], argv[4]);
  printf("Done. Let's eval the circuit...\n");
  executor.Eval(argv[5]);
  printf("Done! Decrypt the result via client.\n");
}
