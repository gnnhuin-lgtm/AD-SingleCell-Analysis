test_that("normalize_data works with CPM method", {
  # 创建测试数据
  test_counts <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 2, ncol = 3)
  
  # 测试CPM标准化
  result <- normalize_data(test_counts, method = "cpm")
  
  # 验证结果
  expect_is(result, "matrix")
  expect_equal(dim(result), dim(test_counts))
  expect_true(all(result >= 0))
})

test_that("log_transform works correctly", {
  test_data <- matrix(c(1, 10, 100), nrow = 3, ncol = 1)
  result <- log_transform(test_data, pseudocount = 1)
  
  expect_is(result, "matrix")
  expect_equal(dim(result), dim(test_data))
})
