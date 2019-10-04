/*
    g++ client_encrypt.cc -o client_encrypt -O3 -I ./tfhe/src/include/ -L ./tfhe/build/libtfhe/ -ltfhe-spqlios-fma
 */

#include <cstdio>
#include <cassert>
#include <string>
#include <array>
#include <tfhe.h>
#include <tfhe_io.h>

#include "constants.h"

using namespace std;

int main(int argc, char *argv[]) {
  if (argc < 5) {
    assert(argc >= 1);
    printf("Usage: %s (query_path) (path_for_secret_key) (path_for_cloud_key) (path_for_enc_query)\n", argv[0]);
    exit(0);
  }

  puts("Setting parameter...");
  TFheGateBootstrappingParameterSet* params = new_default_gate_bootstrapping_parameters(minimum_lambda);

  puts("Setting seed...");
  FILE *rand_dev = fopen("/dev/urandom", "rb");
  array<uint32_t, 3> seed;
  if (!rand_dev) {
    puts("Warning: /dev/urandom is not found. We will use the default seed.");
    seed = { 314, 1592, 657 };
  } else {
    fread(&seed[0], sizeof(uint32_t), 3, rand_dev);
  }
  tfhe_random_generator_setSeed(&seed[0], 3);

  puts("Generating key...");
  TFheGateBootstrappingSecretKeySet* key = new_random_gate_bootstrapping_secret_keyset(params);

  FILE* secret_key = fopen(argv[2],"wb");
  if (!secret_key) {
    printf("Could not create the file %s.\n", argv[1]);
    exit(-1);
  }
  export_tfheGateBootstrappingSecretKeySet_toFile(secret_key, key);
  fclose(secret_key);

  FILE* cloud_key = fopen(argv[3],"wb");
  if (!cloud_key) {
    printf("Could not create the file %s.\n", argv[3]);
    exit(-1);
  }
  export_tfheGateBootstrappingCloudKeySet_toFile(cloud_key, &key->cloud);
  fclose(cloud_key);

  puts("Opening query...");
  FILE *query_fp = fopen(argv[1], "rb");
  if (!query_fp) {
    printf("Could not open the file %s.\n", argv[1]);
    exit(-1);
  }

  string query = "";
  bool suspecting_corrupted = false;
  while (1) {
    int d = fgetc(query_fp);
    if (d == EOF) break;
    
    char c = d;
    if (c == '0' || c == '1') query += c;
    else if (c != '\t' && !suspecting_corrupted) {
      suspecting_corrupted = true;
      printf("Warning: the query file %s seems corrupted.\n", argv[1]);
    }
  }
  assert(query.size() == LENGTH);

  puts("Encrypting query...");
  LweSample* ct = new_gate_bootstrapping_ciphertext_array((int)query.size(), params);
  for (int i=0; i<(int)query.size(); i++) {
    bootsSymEncrypt(&ct[i], query[i] == '1', key);
  }

  puts("Saving encrypted query...");
  FILE* cloud_data = fopen(argv[4],"wb");
  if (!cloud_data) {
    printf("Could not create the file %s.\n", argv[4]);
    exit(-1);
  }
  for (int i=0; i<(int)query.size(); i++) {
    export_gate_bootstrapping_ciphertext_toFile(cloud_data, &ct[i], params);
  }
  fclose(cloud_data);

  delete_gate_bootstrapping_ciphertext_array((int)query.size(), ct);
  delete_gate_bootstrapping_secret_keyset(key);
  delete_gate_bootstrapping_parameters(params);
}
