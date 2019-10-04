# NISMM
This is the experimental implementation of the method of finding set-maximal matches over non-interactive secure computation

# Requirement

 - TFHE
 - sdsl-lite

# How to use
  
  1. Generate keys and encrypt a query with client_encrypt.cc
  2. Preprocess a dictionary with preprocess.cc (it parses vcf)
  3. Run server.cc
  4. Decrypt the result with client_decrypt.cc

  Please modify THRESHOLD, LENGTH, minimum_lambda in constants.h before using those program. 
