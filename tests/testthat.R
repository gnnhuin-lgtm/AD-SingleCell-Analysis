# 该文件被testthat::test_local()调用
if (requireNamespace("testthat", quietly = TRUE)) {
  testthat::test_check("bioinformaticsSkills")
}