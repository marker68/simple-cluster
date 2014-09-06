/*
 * test_kmeans.cpp
 *
 *  Created on: 2014/09/06
 *      Author: Nguyen Anh Tuan <t_nguyen@hal.t.u-tokyo.ac.jp>
 */

#include <iostream>
#include <vector>
#include <gtest/gtest.h>
#include "k-means.h"

using namespace std;
using namespace SimpleCluster;

/**
 * Customized test case for testing
 */
class KmeansTest : public ::testing::Test {
protected:
  // Per-test-case set-up.
  // Called before the first test in this test case.
  // Can be omitted if not needed.
  static void SetUpTestCase() {
    size = (int *)malloc(10 * sizeof(int));
    size[0] = 100;
    for(int i=1;i<10;i++) {
	size[i] = 10 * size[i-1];
    }
  }

  // Per-test-case tear-down.
  // Called after the last test in this test case.
  // Can be omitted if not needed.
  static void TearDownTestCase() {
    delete size;
    size = NULL;
  }

  // You can define per-test set-up and tear-down logic as usual.
  virtual void SetUp() { }
  virtual void TearDown() {}

public:
  // Some expensive resource shared by all tests.
  static int * size;
};

int * KmeansTest::size = NULL;

TEST_F(KmeansTest, test1) {

}

/*int main(int argc, char * argv[]) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}*/
