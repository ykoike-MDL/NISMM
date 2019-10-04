/*
    g++ client_decrypt.cc -o client_decrypt -O3 -I ./tfhe/src/include/ -L ./tfhe/build/libtfhe/ -ltfhe-spqlios-fma
 */

#include <cstdio>
#include <cassert>
#include <string>
#include <array>
#include <algorithm>
#include <tfhe.h>
#include <tfhe_io.h>
#include "constants.h"
using namespace std;


int main(int argc, char *argv[]) {
  setbuf(stdout, NULL);

  if (argc < 4) {
    assert(argc >= 1);
    printf("Usage: %s (secret_key_path) (enc_answer_path) (delimiter_path)\n", argv[0]);
    exit(0);
  }

  FILE* secret_key = fopen(argv[1],"rb");
  if (!secret_key) {
    printf("Could not open the secret key %s.\n", argv[1]);
    exit(-1);
  }

  TFheGateBootstrappingSecretKeySet* key = new_tfheGateBootstrappingSecretKeySet_fromFile(secret_key);
  fclose(secret_key);

  const TFheGateBootstrappingParameterSet* params = key->params;

  puts("Reading the delimiter information...");
  FILE *delim_fp = fopen(argv[3], "rb");
  if (!delim_fp) {
    printf("Could not open the delimiter %s.\n", argv[3]);
    exit(-1);
  }

  auto u16 = [&delim_fp, &argv]() {
    uint16_t val;
    int res = fread(&val, sizeof(uint16_t), 1, delim_fp);
    if (res != 1) {
      printf("The delimiter file %s is corrupted\n", argv[3]);
      exit(-1);
    }
    return val;
  };

  auto u32 = [&delim_fp, &argv]() {
    uint32_t val;
    int res = fread(&val, sizeof(uint32_t), 1, delim_fp);
    if (res != 1) {
      printf("The delimiter file %s is corrupted\n", argv[3]);
      exit(-1);
    }
    return val;
  };

  auto len = u32();
  LweSample* answer = new_gate_bootstrapping_ciphertext_array(len, params);

  auto bitL = u16();
  auto N = u32();
  fclose(delim_fp);

  puts("Done. Reading the encrypted answer...");
  FILE* answer_data = fopen(argv[2], "rb");
  if (!answer_data) {
    printf("Could not open the answer %s.\n", argv[2]);
    exit(-1);
  }

  for (auto i=0; i<len; i++) {
    import_gate_bootstrapping_ciphertext_fromFile(answer_data, &answer[i], params);
  }
  fclose(answer_data);

  puts("Done. Decrypting the answer...");

  uint32_t idx = 0;
  using Seg = pair<int, int>;
  vector<Seg> ls;
  
  for (auto i=0; i<=N; i++) {
    int L = 0;
    for (auto j=0; j<bitL; j++) {
      int ai = bootsSymDecrypt(&answer[idx++], key);
      L |= ai << j;
    }
    if (L == 0) continue;
    ls.emplace_back(Seg(LENGTH-1-THRESHOLD-i, LENGTH-1-THRESHOLD-i+L));
  }
  sort(ls.begin(), ls.end());

  printf("The number of SMM: %d\n", (int)ls.size());
  printf("The lengths:\n");
  for (auto &p : ls) {
    printf("[%d, %d)\n", p.first, p.second);
  }

  delete_gate_bootstrapping_ciphertext_array(len, answer);
  delete_gate_bootstrapping_secret_keyset(key);
}
