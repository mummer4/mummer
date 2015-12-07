// Copyright 2006, Google Inc.
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
//     * Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above
// copyright notice, this list of conditions and the following disclaimer
// in the documentation and/or other materials provided with the
// distribution.
//     * Neither the name of Google Inc. nor the names of its
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#include <iostream>
#include <chrono>
#include <cstring>


#include "gtest/gtest.h"
#include "gtest/test.hpp"

auto seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
std::mt19937 rand_gen;

void get_seed(int& argc, char** argv) {
  int i = 1;
  for( ; i < argc; ++i) {
    if(!strcmp(argv[i], "--seed")) {
      ++i;
      break;
    }
  }
  if(i == argc) return;
  seed = std::atoll(argv[i]);
  for(++i; i < argc; ++i)
    argv[i - 2] = argv[i];
  argc -= 2;
}

std::string sequence(size_t len) {
  static char bases[4] = { 'a', 'c', 'g', 't' };
  std::uniform_int_distribution<int> dist(0, 3);
  std::string                        result;
  for(size_t i = 0; i < len; ++i)
    result += bases[dist(rand_gen)];
  return result;
}

GTEST_API_ int main(int argc, char **argv) {
  get_seed(argc, argv);
  rand_gen.seed(seed);
  std::cout << "Seed:" << seed << std::endl;
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
