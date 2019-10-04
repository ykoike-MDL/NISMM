/*
    g++ preprocess.cc -o preprocess -O3 -I ./sdsl-lite/include/ -L ./sdsl-lite/lib -lsdsl -ldivsufsort -ldivsufsort64
 */

#include <cstdio>
#include <cassert>
#include <iostream>
#include <fstream>
#include <vector>
#include <set>
#include <array>
#include <unordered_map>

#include <sdsl/suffix_arrays.hpp>
#include <sdsl/suffix_trees.hpp>

#include "constants.h"

using namespace std;
using namespace sdsl;

using Size = uint32_t;
using VC = vector<char>;

const char g_num_fanin[18] = {-1, 0, 0, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3};

pair<vector<VC>, VC> ParseTestcase(char *path) {
  char buf[114514];
  FILE *fp = fopen(path, "rb");
  assert (fp);

  int N, M;
  VC query;
  vector<VC> ret;

  fscanf(fp, "%d%d", &N, &M);
  for (int i=0; i<=M; i++) {
    fscanf(fp, " %s", buf);
    assert(strlen(buf) == N);
    VC vc(N);
    for (int j=0; j<N; j++) {
      vc[j] = buf[j] - '0';
    }
    if (i < M) ret.emplace_back(vc);
    else query.swap(vc);
  }
  return make_pair(ret, query);
}

// Simple Parser for VCF. May be inaccurate.
vector<VC> ParseVcf(char *vcf_path) {
  ifstream ifs(vcf_path);
  if (ifs.fail()) {
    printf("Could not open the vcf file %s.\n", vcf_path);
    exit(-1);
  }

  string buf;
  // Skip lines starting with '##'
  while (1) {
    if (!getline(ifs, buf)) {
      printf("Could not parse the vcf file %s.\n", vcf_path);
      exit(-1);
    }

    if (buf.size() < 2) break;
    if (buf[0] != '#' || buf[1] != '#') break;
  }

  // Parse a header line
  if (buf.empty() || buf[0] != '#') {
    printf("Could not parse the vcf file %s.\n", vcf_path);
    printf("It doesn't contain the header line.\n");
    exit(-1);
  }

  size_t idx = 1;
  size_t num_skip = 0;
  set<string> known_column{
    "CHROM",
    "POS",
    "ID",
    "REF",
    "ALT",
    "QUAL",
    "FILTER",
    "INFO",
    "FORMAT"
  };

  // Skip unnecessary columns
  while (1) {
    auto pos = buf.find("\t", idx);
    if (pos == string::npos) pos = buf.size();
    auto pos2 = buf.substr(idx, pos-idx).find(" ", 0); // seems sometimes the file uses ' ' instead of '\t' ...
    if (pos2 != string::npos) {
      assert(idx + pos2 < pos);
      pos = idx + pos2;
    }

    string col = buf.substr(idx, pos-idx);

    if (col.empty()) { // ' ' is used... f***
      idx = pos+1;
      continue;
    }

    if (!known_column.count(col)) {
      if (col != "0") {
        printf("Warning: the column %s is not expected to appear.\n", col.c_str());
        puts("ParseVcf will probably fail parsing the file.");
      }
      break;
    }

    if (pos >= buf.size()) {
      printf("Could not parse the vcf file %s.\n", vcf_path);
      printf("The actual strings are not found.\n");
      exit(-1);
    }
    idx = pos+1;
    ++num_skip;
  }

  // Count the number of string
  size_t num_str = 0;
  while (1) {
    auto pos = buf.find("\t", idx);
    ++num_str;
    if (pos == string::npos) break;
    idx = pos+1;
  }

  vector<VC> ret(num_str*2);
  while (getline(ifs, buf)) {
    size_t idx = 0;
    for (auto i=0; i<num_skip; i++) {
      auto pos = buf.find("\t", idx);
      if (pos == string::npos) {
        printf("Could not parse the vcf file %s.\n", vcf_path);
        printf("The number of column is inconsistent with the previous lines.\n");
        exit(-1);
      }
      idx = pos+1;
    }

    for (auto i=0; i<num_str; i++) {
      auto pos = buf.find("\t", idx);
      if (pos == string::npos) pos = buf.size();

      string val = buf.substr(idx, pos-idx);
      if (val.size() < 3 || val[1] != '|') {
        printf("Could not parse the vcf file %s.\n", vcf_path);
        printf("'%s' doesn't have a form like '0|1'\n", val.c_str());
        exit(-1);
      }

      assert(val[0] == '0' || val[0] == '1');
      assert(val[2] == '0' || val[2] == '1');

      ret[i*2].emplace_back(val[0] == '1');
      ret[i*2+1].emplace_back(val[2] == '1');
      
      idx = pos+1;
    }
  }

  return ret;
}

class Automaton {
private:
  using Gate = tuple<LogicType, Size, Size>;

  struct GateHash : public unary_function<Gate, size_t> {
    std::size_t operator()(const Gate& k) const {
      return std::get<0>(k) ^ std::get<1>(k) ^ std::get<2>(k);
    }
  };

  int bitM_;
  int bitL_;
  Size N_;
  Size M_;
  Size K_;
  Size num_logic_;
#ifdef DEBUG
  TypeForM num_last_states_;
  TypeForM breadth_last_report_;
  vector<TypeForM> breadth_report_;
  vector<vector<TypeForM>> last_reported_;
  vector<array<vector<vector<TypeForM>>, 2>> reported_;
#endif 
  vector<array<vector<TypeForM>, 2>> dest_;
  vector<array<vector<Size>, 2>> reset_length_;
  cst_sct3<csa_bitcompressed<int_alphabet<>>> cst_;
  unordered_map<Gate, Size, GateHash, equal_to<Gate>> gate_cache_;
  FILE *delim_fp_;
  FILE *database_fp_;
  FILE *indices_fp_;
  char *delim_path_;
  char *database_path_;
  char *indices_path_;

  vector<Size> CollectSubtreeIds (Size idx) {
    vector<Size> res;
    const auto root = cst_.root();
    for (char b=0; b<2; b++) {
      auto gc = 1 + idx*2 + (Size)b;
      auto subt = cst_.child(root, gc);
      if (subt == root) continue;
      for_each(cst_.begin(subt), cst_.end(subt), 
        [this, &res] (decltype(root) v) {
          res.emplace_back(cst_.id(v));
        }
      );
    }
    sort(res.begin(), res.end());
    res.erase(unique(res.begin(), res.end()), res.end());
    /*cout << "size: " << res.size() << endl;
      for (auto idx : res) {
      cout << idx << " ";
      }cout << endl;*/
    return res;
  };


  inline vector<Size> CreateMultiMux(Size a, const vector<Size> &b, const vector<Size> &c) {
    auto n = b.size();
    assert(b.size() == c.size());
    assert(!b.empty());

    vector<Size> ret(n);
    for (auto i=0; i<n; i++) {
      ret[i] = num_logic_++;
      PackDatabase8((uint8_t)MUX);
      PackDatabase32(a);
      PackDatabase32(b[i]);
      PackDatabase32(c[i]);
    }
    return ret;
  }

  inline Size CreateGate(LogicType type, const Size &a, const Size &b) {
    auto g = Gate{type, a, b};
    auto itr = gate_cache_.find(g);
    if (itr != gate_cache_.end()) return itr->second;

    auto new_idx = num_logic_++;
    gate_cache_[g] = new_idx;
    PackDatabase8((uint8_t)type);
    PackDatabase32(a);
    PackDatabase32(b);
    return new_idx;
  }

  inline Size CreateNotGate(const Size &a) {
    auto g = Gate{NOT, a, a};
    auto itr = gate_cache_.find(g);
    if (itr != gate_cache_.end()) return itr->second;
    
    auto new_idx = num_logic_++;
    gate_cache_[g] = new_idx;
    PackDatabase8((uint8_t)NOT);
    PackDatabase32(a);
    return new_idx;
  }

  inline Size CreateMultiGate(LogicType type, const vector<Size> &vs) {
    auto n = vs.size();
    if (n == 0) {
      return N_;
    }

    if (n == 1) {
      return vs[0];
    }

    using Tup = pair<Size, Size>;
    priority_queue<Tup, vector<Tup>, greater<Tup>> q;

    for (auto i=0; i<n; i++) {
      q.push(Tup(1, vs[i]));
    }

    while (q.size() > 1) {
      auto a = q.top(); q.pop();
      auto b = q.top(); q.pop();
      Size new_idx;

      new_idx = CreateGate(type, a.second, b.second);
      q.push(Tup(a.first+b.first, new_idx));
    }
    assert(q.size() == 1);
    return q.top().second;
  }

  Size CreateClause(const vector<Size> &vs, const vector<char> &neg) {
    assert(neg.size() == vs.size());
    auto n = vs.size();
    assert(n >= 1);
    if (n == 1) {
      if (neg[0]) return CreateNotGate(vs[0]);
      return vs[0];
    }

    using Tup = tuple<Size, Size, char>;
    priority_queue<Tup, vector<Tup>, greater<Tup>> q;

    for (auto i=0; i<n; i++) {
      q.push(Tup{1, vs[i], neg[i]});
    }

    while (q.size() > 1) {
      auto a = q.top(); q.pop();
      auto b = q.top(); q.pop();
      Size new_idx;
      if (get<2>(a) && get<2>(b)) {
        new_idx = CreateGate(NOR, get<1>(a), get<1>(b));
      } else if (get<2>(a)) {
        new_idx = CreateGate(ANDNY, get<1>(a), get<1>(b));
      } else if (get<2>(b)) {
        new_idx = CreateGate(ANDYN, get<1>(a), get<1>(b));
      } else {
        new_idx = CreateGate(AND, get<1>(a), get<1>(b));
      }
      q.push(Tup{get<0>(a)+get<0>(b), new_idx, 0});
    }
    assert(q.size() == 1);
    return get<1>(q.top());
  }

  vector<char> BitResolve(Size val, const size_t width=0) {
    vector<char> ret;
    while (val > 0) {
      ret.emplace_back(val&1);
      val >>= 1;
    }
    
    if (width) {
      assert(ret.size() <= width);
      ret.resize(width);
    }
    return ret;
  }

  int CalcBitSize(Size t) {
    int ret = 0;
    for (Size k=1; k<t; k <<= 1) ++ret;
    return max(1, ret);
  }

  vector<Size> CreateMinterms(const set<int> needed, const vector<Size> &bits) {
    if (bits.empty()) {
      if (!needed.empty()) {
        assert(*needed.begin() == 0);
      }
      return vector<Size>(1, N_+1);
    }

    int n = bits.size();
    assert(0 < n && n <= 30);
    map<long long, Size> mp;
    vector<Size> ret(1 << n);

    for (int i=0; i<n; i++) {
      mp[1LL << (i*2)] = bits[i];
    }

    const Size INF = -1;
    assert(num_logic_ <= INF);
    function<Size(int, int, int)> Get;
    Get = [&](long long st, int l, int r) {
      assert(r - l > 0);
      st &= ((1LL << r) - 1);
      st &= ~((1LL << l) - 1);
      if (st == 0) return INF;
      if (mp.count(st)) return mp[st]; 
      if (__builtin_popcountll(st) == 1) {
        int idx = __builtin_ffsll(st);
        assert(st == (1 << (idx-1)));
        assert(idx > 0);
        assert(idx%2 == 0);
        assert((idx-1)/2 >= 0);
        assert(n > (idx-1)/2);
        return mp[st] = CreateNotGate(bits[(idx-1)/2]);
      }
     
      int mid = (l+r)/2;
      Size lv = Get(st, l, mid);
      Size rv = Get(st, mid, r);
      if (lv == INF) {
        assert(mp.count(st));
        assert(mp[st] == rv);
        return rv;
      }
      if (rv == INF) {
        assert(mp.count(st));
        assert(mp[st] == lv);
        return lv;
      }
      return mp[st] = CreateGate(AND, lv, rv);
    };

    for (int st : needed) {
      long long v = 0;
      for (int i=0; i<n; i++) {
        if (st >> i & 1) v |= 1LL << (i*2);
        else v |= 1LL << (i*2+1);
      }
      assert(v > 0);
      ret[st] = Get(v, 0, n*2);
      assert(ret[st] != 0);
    }
    return ret;
  }

  Size CreateMapping(const vector<Size> &bits, int width, const vector<char> &onof) {
    assert(width >= 2);
    assert((int)bits.size() == width);
    assert((int)onof.size() == 1 << width);

    // (k, s)-lupanov representation
    // https://www.math.ucsd.edu/~sbuss/CourseWeb/Math262A_2013F/Scribe02.pdf
    int k = min((int)ceil(3 * log2(width)), width);
    int s = max((int)ceil(width - 5 * log2(width)), 1);
    int num_rows = 1 << k;
    int p = num_rows / s + (num_rows%s != 0);
    int num_cols = 1 << (width-k);
    assert(s < 30);
    int lim = 1 << s;

    set<int> col_minterms;
    set<int> row_minterms;
    vector<vector<vector<int>>> true_rows;
    vector<vector<vector<int>>> true_cols;

    vector<Size> col_bits(width-k);
    vector<Size> row_bits(k);
    for (int i=0; i<width-k; i++) {
      col_bits[i] = bits[i];
    }

    for (int i=0; i<k; i++) {
      row_bits[i] = bits[i+width-k];
    }

    for (int i=0; i<p; i++) { // fix A_i

      // process f_col(b)
      true_rows.emplace_back(vector<vector<int>>(lim, vector<int>()));
      true_cols.emplace_back(vector<vector<int>>(lim, vector<int>()));
      
      for (int j=0; j<num_cols; j++) {
        int v = i*s;
        int b = 0;
        for (int l=0; l<s; l++) {
          if (v >= num_rows) break;
          assert(v * num_cols + j < onof.size());
          b = (b << 1) | onof[v * num_cols + j];
          ++v;
        }
        true_cols[i][b].emplace_back(j);
        col_minterms.insert(j);
      }

      // process f_row
      for (int b=0; b<lim; b++) { // fix b
        int v = i*s;

        for (int j=0; j<s; j++) {
          if (v >= num_rows) break;
          if (b >> j & 1) {
            true_rows[i][b].emplace_back(v);
            row_minterms.insert(v);
          }
          ++v;
        } 
      }
    }

    vector<Size> clauses;
    auto col_map = CreateMinterms(col_minterms, col_bits);
    auto row_map = CreateMinterms(row_minterms, row_bits);
    for (int i=0; i<p; i++) { // fix A_i
      for (int b=0; b<lim; b++) {
        if (true_cols[i][b].empty()) continue;
        if (true_rows[i][b].empty()) continue;
        
        vector<Size> col_dnf_cls;
        for (int v : true_cols[i][b]) {
          assert(col_map.size() > v);
          assert(col_map[v] != 0);
          col_dnf_cls.emplace_back(col_map[v]);
        }

        vector<Size> row_dnf_cls;
        for (int v : true_rows[i][b]) {
          assert(row_map[v] != 0);
          row_dnf_cls.emplace_back(row_map[v]);
        }
        
        Size cold = CreateMultiGate(OR, col_dnf_cls);
        Size rowd = CreateMultiGate(OR, row_dnf_cls);
        clauses.emplace_back(CreateGate(AND, cold, rowd));
      }
    }
    return CreateMultiGate(OR, clauses);
  }

public:
  void DebugMapping(const char *database_path, const char *indices_path) {
    database_path_ = (char*)database_path;
    database_fp_ = fopen(database_path, "wb");
    if (!database_fp_) {
      printf("Error: couldn't open the database file %s.\n", database_path);
      exit(-1);
    }

    indices_path_ = (char*)indices_path;
    indices_fp_ = fopen(indices_path, "wb");
    if (!indices_fp_) {
      printf("Error: couldn't open the indices file %s.\n", indices_path);
      exit(-1);
    }

    K_ = 0;

    N_ = 5;

    PackDatabase32(0xdeadbeefU);
    PackDatabase32(N_);

    vector<Size> bits(N_);
    for (int i=0; i<N_; i++) {
      PackDatabase8(INPUT);
      bits[i] = i;
    }
    PackDatabase8(CONSTANT);
    PackDatabase8(CONSTANT);
    num_logic_ = N_ + 2;

    vector<char> onof(1 << N_, 0);
    for (int i=0; i<(1 << N_); i++) {
      onof[i] = __builtin_popcountll(i)%2;
    }
    Size out = CreateMapping(bits, N_, onof);
    PackIndices32(out);
    EndIndices();
    PackDatabase8(ERROR);

    fclose(database_fp_);
    fclose(indices_fp_);
  }

  vector<Size> CreateTransition(Size idx, const vector<Size> &state_bits) {
    vector<Size> bits(1, idx-1);
    bits.insert(bits.end(), state_bits.begin(), state_bits.end());
    vector<char> onof;

    auto num_states = dest_[idx][0].size();
    int A = CalcBitSize(num_states) + 1;
    assert(A <= 16);
    assert((int)state_bits.size()+1 == A);

    Size num_nxt_states = idx > 1 ? dest_[idx-1][0].size() : num_last_states_;
    int B = CalcBitSize(num_nxt_states);
    vector<Size> nxt_state_bits(B);
    for (int r=0; r<B; r++) {
      onof.assign(1 << bits.size(), 0);

      for (char b=0; b<2; b++) {
        assert(dest_[idx][b].size() == num_states);
        for (auto i=0; i<num_states; i++) {
          if (dest_[idx][b][i] >> r & 1) {
            onof[i*2 + (Size)b] = 1;
          } 
        }
      }

      nxt_state_bits[r] = CreateMapping(bits, A, onof);
    }

    return nxt_state_bits;
  }

  Size CreateResetJudge(Size idx, const vector<Size> &state_bits) {
    vector<Size> bits(1, idx-1);
    bits.insert(bits.end(), state_bits.begin(), state_bits.end());

    vector<char> onof;

    auto num_states = dest_[idx][0].size();
    int A = CalcBitSize(num_states) + 1;
    assert(A <= 16);
    assert((int)state_bits.size()+1 == A);

    onof.assign(1 << bits.size(), 0);
    for (char b=0; b<2; b++) {
      assert(dest_[idx][b].size() == num_states);
      assert(reset_length_[idx][b].size() == num_states);
      for (auto i=0; i<num_states; i++) {
        if (reset_length_[idx][b][i] != N_+1) {
          onof[i*2 + (Size)b] = 1;
        } 
      }
    }

    return CreateMapping(bits, A, onof);
  }

  vector<Size> CreateReset(Size idx, const vector<Size> &state_bits, int bitL) {
    vector<Size> bits(1, idx-1);
    bits.insert(bits.end(), state_bits.begin(), state_bits.end());

    vector<char> onof;

    auto num_states = dest_[idx][0].size();
    int A = CalcBitSize(num_states) + 1;
    assert(A <= 16);
    assert((int)state_bits.size()+1 == A);

    vector<Size> ret(bitL);
    for (int r=0; r<bitL; r++) {
      onof.assign(1 << bits.size(), 0);

      for (char b=0; b<2; b++) {
        assert(reset_length_[idx][b].size() == num_states);
        for (auto i=0; i<num_states; i++) {
          auto l = reset_length_[idx][b][i];
          auto v = i*2 + (Size)b;
          if (l != N_+1 && (l >> r & 1)) {
            onof[v] = 1;
          }
        }
      } 

      ret[r] = CreateMapping(bits, A, onof);
    }

    return ret;
  }

  // FIXME: can be faster
  vector<Size> CreateIncrementer(Size num, const vector<Size> &L_bits) {
    int n = L_bits.size();
    assert(n == bitL_);
    vector<Size> ret(n);
    ret[0] = CreateNotGate(L_bits[0]);
    vector<Size> lowers(1, L_bits[0]);
    lowers.reserve(n);
    for (int i=1; i<n; i++) {
      ret[i] = CreateGate(XOR, CreateMultiGate(AND, lowers), L_bits[i]);
      if (i < n-1) lowers.emplace_back(L_bits[i]);
    }
    return ret;
  }

  // FIXME: can be faster
  Size CreateThresholdJudge(const vector<Size> &L_bits) {
    int lim = CalcBitSize(THRESHOLD);
    int m = (int)L_bits.size() - lim;
    assert((int)L_bits.size() == bitL_);
    assert(m >= 0);

    vector<Size> mor_inputs(m);
    for (int i=0; i<m; i++) {
      mor_inputs[i] = L_bits[lim+i];
    }
    mor_inputs.reserve(m + lim);

    vector<Size> is_eqs;
    for (int i=lim-1; i>=0; i--) {
      if (THRESHOLD >> i & 1) {
        is_eqs.emplace_back(L_bits[i]);
      } else {
        if (is_eqs.empty()) {
          mor_inputs.emplace_back(L_bits[i]);
        } else if (i > 0) {
          auto cond = CreateMultiGate(AND, is_eqs);
          mor_inputs.emplace_back(CreateGate(AND, cond, L_bits[i]));
        }
      }
    }

    if (!is_eqs.empty()) {
      return CreateGate(OR, CreateMultiGate(OR, mor_inputs), CreateMultiGate(AND, is_eqs));
    }

    return CreateMultiGate(OR, mor_inputs);
  }

  void CreateOutput(Size idx, Size should_report, const vector<Size> &L_bits, const vector<Size> &state_bits) {
    int n = (int)L_bits.size();
    assert(n == bitL_);
    for (int i=0; i<n; i++) {
      PackIndices32(CreateGate(AND, should_report, L_bits[i]));
    }
  }

  void PackDatabase8(uint8_t val) {
    auto res = fwrite(&val, sizeof(uint8_t), 1, database_fp_);
    if (res != 1) {
      printf("Could not write a value in %s...\n", database_path_);
      exit(-1);
    }
  }

  void PackDatabase16(uint16_t val) {
    auto res = fwrite(&val, sizeof(uint16_t), 1, database_fp_);
    if (res != 1) {
      printf("Could not write a value in %s...\n", database_path_);
      exit(-1);
    }
  }

  void PackDatabase32(uint32_t val) {
    auto res = fwrite(&val, sizeof(uint32_t), 1, database_fp_);
    if (res != 1) {
      printf("Could not write a value in %s...\n", database_path_);
      exit(-1);
    }
  }

  void PackIndices32(uint32_t val) {
    ++val;
    ++K_;
    auto res = fwrite(&val, sizeof(uint32_t), 1, indices_fp_);
    if (res != 1) {
      printf("Could not write a value in %s...\n", indices_path_);
      exit(-1);
    }
  }

  void EndIndices() {
    uint32_t val = 0;
    auto res = fwrite(&val, sizeof(uint32_t), 1, indices_fp_);
    if (res != 1) {
      printf("Could not write a value in %s...\n", indices_path_);
      exit(-1);
    }
  }

  void PackDelim8(uint8_t val) {
    auto res = fwrite(&val, sizeof(uint8_t), 1, delim_fp_);
    if (res != 1) {
      printf("Could not write a value in %s...\n", delim_path_);
      exit(-1);
    }
  }

  void PackDelim16(uint16_t val) {
    auto res = fwrite(&val, sizeof(uint16_t), 1, delim_fp_);
    if (res != 1) {
      printf("Could not write a value in %s...\n", delim_path_);
      exit(-1);
    }
  }

  void PackDelim32(uint32_t val) {
    auto res = fwrite(&val, sizeof(uint32_t), 1, delim_fp_);
    if (res != 1) {
      printf("Could not write a value in %s...\n", delim_path_);
      exit(-1);
    }
  }

public:
  Automaton() {}

  void Init(const vector<VC> &str_set) {
    M_  = str_set.size();
    if (M_ < 1 || 1000 < M_) {
      puts("Error: currently we consider only up to M=1000.");
      exit(-1);
    }
    bitM_ = CalcBitSize(M_);
    assert(bitM_ == CalcBitSize(M_+1));

    if (M_*2 > (1ULL << LIM_BITM)) {
      puts("Error: currently we consider only up to M=1000.");
      puts("... Change TypeForM to uint32_t...");
      exit(-1);
    }

    N_ = str_set[0].size();
    assert(N_ == LENGTH);
    if (N_*M_ < N_) {
      puts("Error: no way! You gave a really too large dataset...");
      puts("We cannot guarantee this code works correctly in such a large testcase.");
      puts("... Change Size to uint64_t...");
      exit(-1);
    }

    printf("  Done. Constructing CST... This may take 5 minutes or so...\n");
    //if (!load_from_file(cst_, cache_path)) {
    int bitN = 0;
    Size sigma = 1 + N_*2 + M_;
    for (Size k=1; k<sigma; k <<= 1) bitN++;
    printf("  Joining all of the strings...\n");
    Size joined_sz = N_*M_ + M_;
    int_vector<> joined_str(joined_sz, 0, bitN);
    Size idx = 0;
    for (auto j=0; j<M_; j++) {
      assert(str_set.size() > j);
      auto &str = str_set[j];
      assert(str.size() == N_);
      for (auto i=0; i<N_; i++) {
        assert(str[i] == 0 || str[i] == 1);
        joined_str[idx++] = 1 + i*2 + (Size)str[i];
      }
      joined_str[idx++] = 1 + N_*2 + j;
    }
    assert(idx == joined_sz);

    construct_im(cst_, joined_str);
    //store_to_file(cst_, cache_path);
    //}
    printf("  Done. Now, let's create an automaton! This may take 45 minutes or so...\n");

    const auto root = cst_.root();
    using Node = typename remove_const<decltype(root)>::type;

    dest_.resize(N_);
    reset_length_.resize(N_);
#ifdef DEBUG
    breadth_report_.resize(N_);
    reported_.resize(N_);
#endif 

    auto dict = CollectSubtreeIds(N_-1);
    if (dict.size()+1 >= (1ULL << LIM_BITM)) {
      puts("Error: no way, strings are too many! Please change TypeForM to uint32_t....");
      exit(-1);
    }

    auto ReverseDict = [this] (const decltype(dict) &dic, const Node &v) {
      auto idx = cst_.id(v);
      auto itr = lower_bound(dic.begin(), dic.end(), idx);
      assert(*itr == idx);
      return itr - dic.begin();
    };

    for (char b=0; b<2; b++) {
      dest_[0][b].resize(1);
      reset_length_[0][b].resize(1);
      Size qc = 1 + (N_-1)*2 + (Size)b;
      auto subt = cst_.child(root, qc);
      if (subt == root) {
        assert(dict.size() + 1 == dest_[N_-1][b].size());
        dest_[0][b][0] = dict.size();
        reset_length_[0][b][0] = 0;
      } else {
        auto nidx = ReverseDict(dict, subt);
        dest_[0][b][0] = nidx;
        reset_length_[0][b][0] = N_ + 1;
      }
    }

    auto span = N_/100.0;
    //cout << "N_:" << N_ << endl;
    for (Size idx=N_-1; idx>0; --idx) {
      if ((N_-idx)%100 == 0) printf(" %.3f%%...\n", (N_-idx)/span);

      auto nxt_dict = CollectSubtreeIds(idx-1);
      //cout << "nxt:" << nxt_dict.size() << endl;
      if (nxt_dict.size()+1 >= (1ULL << LIM_BITM)) {
        puts("Error: no way, strings are too many! Please change TypeForM to uint32_t....");
        exit(-1);
      }

      for (char b=0; b<2; b++) {
        dest_[idx][b].resize(dict.size()+1);
        reset_length_[idx][b].resize(dict.size()+1);
#ifdef DEBUG
        reported_[idx][b].resize(dict.size()+1);
#endif
      }

      for (auto i=0; i<dict.size(); i++) {
        const auto vid = dict[i];
        const auto v = cst_.inv_id(vid);
        //cout << "vid: " << vid << endl;
        //cout << "v's id: " << cst_.id(v) << endl;
        for (char b=0; b<2; b++) {
          Size qc = 1 + (idx-1)*2 + (Size)b;
          auto nv = cst_.wl(v, qc);
          if (nv != root) {
            //cout << v << " + " << qc << " = " << nv << "(no reset)" << endl;
            auto nidx = ReverseDict(nxt_dict, nv);
            dest_[idx][b][i] = nidx;
            reset_length_[idx][b][i] = N_+1;
          } else {
            auto iv = root; // this is a sentinel
            auto nump = cst_.csa.bwt.rank(v.i, qc);
            if (nump > 0) {
              auto i = cst_.csa.bwt.select(nump, qc);
              iv = cst_.select_leaf(i+1);
            }
            auto jv = root; // likewise
            auto numq = cst_.csa.bwt.rank(v.j+1, qc);
            auto j = cst_.csa.bwt.select(numq+1, qc);
            if (j < cst_.size()) { // found
              jv = cst_.select_leaf(j+1);
            }
            auto u = cst_.lca(v, iv);
            auto w = cst_.lca(v, jv);
            if (cst_.depth(u) < cst_.depth(w)) u = w;

            auto nv = cst_.wl(u, qc);
            if (nv == root) {
              //cout << v << " + " << qc << " = root(reset to 0)" << endl;
              dest_[idx][b][i] = nxt_dict.size();
              reset_length_[idx][b][i] = 0;
            } else {
              //cout << v << " + " << qc << " = " << nv << "(reset to " << cst_.depth(u) + 1 << ")" <<  endl;
              auto nidx = ReverseDict(nxt_dict, nv);
              dest_[idx][b][i] = nidx;
              reset_length_[idx][b][i] = cst_.depth(u) + 1;
            }

#ifdef DEBUG
            if (cst_.depth(v) >= THRESHOLD) {
              TypeForM wid = v.j - v.i + 1;
              reported_[idx][b][i].resize(v.j - v.i + 1);
              breadth_report_[idx] = max<TypeForM>(breadth_report_[idx], wid);
              for (auto x=v.i; x<=v.j; x++) {
                reported_[idx][b][i][x-v.i] = cst_.csa[x]/(N_+1);
              }
            }
#endif
          }
        }
      }

      for (char b=0; b<2; b++) { 
        assert(dest_[idx][b].size() > dict.size());
        assert(reset_length_[idx][b].size() > dict.size());
        Size qc = 1 + (idx-1)*2 + (Size)b;
        auto subt = cst_.child(root, qc);
        if (subt == root) {
          dest_[idx][b][dict.size()] = nxt_dict.size();
          reset_length_[idx][b][dict.size()] = 0;
        } else {
          auto nidx = ReverseDict(nxt_dict, subt);
          dest_[idx][b][dict.size()] = nidx;
          reset_length_[idx][b][dict.size()] = 1;
        }
      }
      dict.swap(nxt_dict);
    }
    assert(dest_.size() == N_);

#ifdef DEBUG
    num_last_states_ = dict.size()+1;
    last_reported_.resize(num_last_states_);
    for (auto i=0; i+1<num_last_states_; i++) {
      auto vid = dict[i];
      auto v = cst_.inv_id(vid);
      if (cst_.depth(v) >= THRESHOLD) {
        last_reported_[i].resize(v.j - v.i + 1);
        TypeForM wid = v.j - v.i + 1;
        breadth_last_report_ = max<TypeForM>(breadth_last_report_, wid);
        for (auto x=v.i; x<=v.j; x++) {
          last_reported_[i][x-v.i] = cst_.csa[x]/(N_+1);
        }
      }
    }
#endif
  }

#ifdef DEBUG
  void Test(char *test_path) {
    puts("Testing...\nOpening query...");
    FILE *query_fp = fopen(test_path, "rb");
    if (!query_fp) {
      printf("Could not open the file %s.\n", test_path);
      exit(-1);
    }

    VC query;
    bool suspecting_corrupted = false;
    while (1) {
      int d = fgetc(query_fp);
      if (d == EOF) break;

      char c = d;
      if (c == '0' || c == '1') query.emplace_back(c-'0');
      else if (c != '\t' && !suspecting_corrupted) {
        suspecting_corrupted = true;
        printf("Warning: the query file %s seems corrupted.\n", test_path);
      }
    }
    Test(query);
  }

  void Test(const VC &query) {
    assert(query.size() == N_);

    Size v = dest_[0][query[N_ - 1]][0];
    Size L = reset_length_[0][query[N_ - 1]][0] == N_ + 1 ? 1 : 0;

    Size cnt = 0;
    for (Size i=N_-1; i>0; i--) {
      //if (N_ - i > THRESHOLD) printf("%u: %u\n", cnt++, L);

      assert(reset_length_[i][query[i-1]].size() > v);
      if (reset_length_[i][query[i-1]][v] != N_+1) {
        if (L >= THRESHOLD) {
          for (auto idx : reported_[i][query[i-1]][v]) {
            printf("%hu %u %u %u\n", idx, i, i+L, L);
          }
        }

        L = reset_length_[i][query[i-1]][v];
        //printf("i: %u v: %u b: %d reseted to %u\n", i, v, query[i-1], L);
      } else {
        ++L;
        //printf("i: %u v: %u b: %d increed to %u\n", i, v, query[i-1], L);
      }

      assert(dest_[i][query[i-1]].size() > v);
      v = dest_[i][query[i-1]][v];
    }

    if (L >= THRESHOLD) {
      assert(last_reported_.size() > v);
      for (auto idx : last_reported_[v]) {
        printf("%hu %u %u %u\n", idx, 0u, 0+L, L);
      }
    }

    puts("Query Done.");
  }
#endif

  bool ExportAsCircuit(const char *database_path, const char *indices_path, const char *delim_path) {
    database_path_ = (char*)database_path;
    database_fp_ = fopen(database_path, "wb");
    if (!database_fp_) {
      printf("Error: couldn't open the database file %s.\n", database_path);
      exit(-1);
    }

    indices_path_ = (char*)indices_path;
    indices_fp_ = fopen(indices_path, "wb");
    if (!indices_fp_) {
      printf("Error: couldn't open the indices file %s.\n", indices_path);
      exit(-1);
    }

    delim_path_ = (char*)delim_path;
    delim_fp_ = fopen(delim_path, "wb");
    if (!delim_fp_) {
      printf("Error: couldn't open the delim file %s.\n", delim_path);
      exit(-1);
    }

    printf("Done. Exporting the database as a circuit...\n");
    printf("This may take about 50 days...\n");
    puts("  Let's retrieve things left...");

    bitM_ = CalcBitSize(M_);
    assert(bitM_ == CalcBitSize(M_+1));

    bitL_ = CalcBitSize(N_+1);

    K_ = 0;

    PackDatabase32(0xdeadbeefU);
    PackDatabase32(N_);

    for (int i=0; i<N_; i++) {
      PackDatabase8(INPUT);
    }
    PackDatabase8(CONSTANT);
    PackDatabase8(CONSTANT);
    num_logic_ = N_ + 2;

    const Size bs[] = {N_, N_+1};
    vector<Size> L_bits(bitL_, bs[0]);
    bool has_not = false;
    Size ini_not;
    if (reset_length_[0][0][0] == N_+1 && reset_length_[0][1][0] == N_+1) {
      L_bits[0] = bs[1];
    } else if (reset_length_[0][0][0] != N_+1 && reset_length_[0][1][0] != N_+1) {
      L_bits[0] = bs[0];
    } else if (reset_length_[0][1][0] == N_+1) {
      L_bits[0] = N_-1;
    } else {
      ini_not = L_bits[0] = CreateNotGate(N_-1);
      has_not = true;
    }

    vector<Size> ini_bits[2];
    vector<Size> state_bits(bitM_);
    for (int k=0; k<bitM_; k++) {
      char a = dest_[0][0][0] >> k & 1;
      char b = dest_[0][1][0] >> k & 1;
      if (a == b) {
        state_bits[k] = bs[a];
      } else if (b) {
        state_bits[k] = N_ - 1;
      } else if (has_not) {
        state_bits[k] = ini_not;
      } else {
        ini_not = state_bits[k] = CreateNotGate(N_-1);
        has_not = true;
      }
    }

    auto span = N_/100.0;
    for (auto idx=N_-1; idx > 0; idx--) {
      printf("%u transitions remaining... current size: %u\n", idx, num_logic_);
      //if ((N_-idx)%100 == 0) printf(" %.3f%%...\n", (N_-idx)/span);

      gate_cache_.clear();
      auto is_reset = CreateResetJudge(idx, state_bits);
      auto mapped = CreateReset(idx, state_bits, bitL_);
      auto inced = CreateIncrementer(N_-idx, L_bits);
      auto over_threshold = CreateThresholdJudge(L_bits);
      auto should_report = CreateGate(AND, is_reset, over_threshold);

      // see below
      if (N_ - idx > THRESHOLD) {
        CreateOutput(idx, should_report, L_bits, state_bits);
      }

      L_bits = CreateMultiMux(is_reset, mapped, inced);
      state_bits = CreateTransition(idx, state_bits);

      for (char b=0; b<2; b++) {
        vector<TypeForM> p;
        vector<vector<TypeForM>> q;
        vector<Size> r;

        dest_[idx][b].swap(p);
        reported_[idx][b].swap(q);
        reset_length_[idx][b].swap(r);
      }

      if (idx == 0) break; // avoid underflow
    }

    // last report
    auto should_report = CreateThresholdJudge(L_bits);
    CreateOutput(N_, should_report, L_bits, state_bits);

    EndIndices();
    PackDatabase8(ERROR);

    PackDelim32((uint32_t)K_);
    PackDelim16((uint16_t)bitL_);

    int32_t num = N_-1-THRESHOLD;
    if (num < 0) num = 0;
    PackDelim32((uint32_t)num);

    fclose(database_fp_);
    fclose(delim_fp_);
    fclose(indices_fp_);
    return true;
  }
};

int main(int argc, char *argv[]) {
#if 0
  Automaton *am = new Automaton();
  am->DebugMapping("database", "indices");
  return 0;
#else
  setbuf(stdout, NULL);
  setbuf(stderr, NULL);

  if (argc < 5) {
    assert(argc >= 1);
    printf("Usage: %s (vcf_path) (path_for_saving_database) (path_for_saving_output_indices) (path_for_saving_delimiter)\n", argv[0]);
    exit(0);
  }

  printf("Parsing the given vcf file...\n");
  auto str_set = ParseVcf(argv[1]);
  //vector<VC> str_set;
  //VC query;
  //tie(str_set, query) = ParseTestcase(argv[1]);
  printf("Done. Initializing Automaton...\n");
  //str_set = {vector<char>{0, 1, 0, 1, 0}, vector<char>{0, 0, 1, 1, 1}, VC{1, 1, 1, 1, 1}, VC{1, 0, 1, 0, 0}};
  
  Automaton *am = new Automaton();
  am->Init(str_set);
  am->Test(query);
  assert(am->ExportAsCircuit(argv[2], argv[3], argv[4]));
  //am->Test("../test/query");
  printf("Done! Execute ./server with the database file.\n");
#endif
}
