# Revised TSCS 2012 reader used by the manuscript analysis.
#
# The defaults reproduce the primary analysis:
#   - full analytic sample;
#   - combined J1/J2 attitude measure;
#   - three ordered attitude categories.

Read_TSCS2012 <- function(
    data_path = "tscs122_v130218_2012性別平等調查外遇的_.csv",
    z_type = c("combine", "j1", "j2"),
    marriage_filter = c("all", "married", "nonmarried"),
    collapse_z = TRUE) {

  z_type <- match.arg(z_type)
  marriage_filter <- match.arg(marriage_filter)

  if (!file.exists(data_path)) {
    stop("找不到資料檔：", data_path)
  }

  tscs.2012.0 <- read.csv(
    data_path,
    header = TRUE,
    sep = ","
  )

  required_variables <- c("a1", "a11", "a10", "j1", "j2", "b1", "k9", "j5")
  missing_variables <- setdiff(required_variables, names(tscs.2012.0))

  if (length(missing_variables) > 0L) {
    stop("資料缺少必要變項：", paste(missing_variables, collapse = ", "))
  }

  tscs.2012.Greenbergdata <- data.frame(
    a1  = tscs.2012.0$a1,
    a11 = tscs.2012.0$a11,
    a10 = tscs.2012.0$a10,
    j1  = tscs.2012.0$j1,
    j2  = tscs.2012.0$j2,
    b1  = tscs.2012.0$b1,
    k9  = tscs.2012.0$k9,
    j5  = tscs.2012.0$j5
  )

  # Complete analytic sample used in the manuscript: n = 1,838.
  tscs.2012.GR <- subset(
    tscs.2012.Greenbergdata,
    k9 %in% 1:2 & j1 %in% 1:4 & j2 %in% 1:4 & j5 %in% 1:4
  )

  # 0 = non-married (single or cohabiting); 1 = ever married.
  tscs.2012.GR$marital_group <- NA_integer_
  tscs.2012.GR$marital_group[tscs.2012.GR$a10 %in% c(1, 2)] <- 0L
  tscs.2012.GR$marital_group[tscs.2012.GR$a10 %in% c(3, 4, 5, 6)] <- 1L

  tscs.2012.GR <- subset(tscs.2012.GR, !is.na(marital_group))

  if (marriage_filter == "married") {
    tscs.2012.GR <- subset(tscs.2012.GR, marital_group == 1L)
  } else if (marriage_filter == "nonmarried") {
    tscs.2012.GR <- subset(tscs.2012.GR, marital_group == 0L)
  }

  # Education categories used by the original program.
  tscs.2012.GR$edu <-
    ifelse(tscs.2012.GR$b1 %in% 1:4,   1L, 0L) +
    ifelse(tscs.2012.GR$b1 %in% 5:9,   2L, 0L) +
    ifelse(tscs.2012.GR$b1 %in% 10:15, 3L, 0L) +
    ifelse(tscs.2012.GR$b1 %in% 16:19, 4L, 0L) +
    ifelse(tscs.2012.GR$b1 %in% 20:21, 5L, 0L)

  tscs.2012.GR$edu1 <- ifelse(
    tscs.2012.GR$edu <= 3L,
    1L,
    ifelse(tscs.2012.GR$edu == 4L, 2L, 3L)
  )

  # Randomized-response answer: 1 = yes, 0 = no.
  Y1.0 <- ifelse(tscs.2012.GR$k9 == 1L, 1L, 0L)

  # Construct the four-category attitude outcome first.
  if (z_type == "j1") {
    ZZ4 <- 5 - tscs.2012.GR$j1
  } else if (z_type == "j2") {
    ZZ4 <- 5 - tscs.2012.GR$j2
  } else {
    z_avg <- (tscs.2012.GR$j1 + tscs.2012.GR$j2) / 2

    # This is the rule that reproduces the manuscript estimates.
    ZZ4 <- ifelse(
      z_avg <= 1.5, 4L,
      ifelse(z_avg == 2.0, 3L,
        ifelse(z_avg <= 3.0, 2L, 1L)
      )
    )
  }

  # Manuscript three-category outcome:
  # Z = 1: situationally acceptable/acceptable
  # Z = 2: unacceptable
  # Z = 3: absolutely unacceptable
  if (collapse_z) {
    ZZ <- ifelse(ZZ4 == 4L, 3L, ifelse(ZZ4 == 3L, 2L, 1L))
  } else {
    ZZ <- ZZ4
  }

  X1 <- ifelse(tscs.2012.GR$a1 == 1L, 1L, 0L)
  X2 <- ifelse(tscs.2012.GR$edu1 == 1L, 0L, 1L)

  chi.1 <- cbind(1, X1, X2)
  colnames(chi.1) <- c("const.", "gender", "education")

  out <- cbind(
    YY = 1,
    Y1.0 = Y1.0,
    ZZ = ZZ,
    chi.1
  )

  cat("樣本數 =", nrow(out), "\n")
  cat("a10 distribution:\n")
  print(table(tscs.2012.GR$a10))
  cat("marital_group distribution (0=nonmarried incl. cohabiting, 1=ever married):\n")
  print(table(tscs.2012.GR$marital_group))
  cat("ZZ distribution:\n")
  print(table(ZZ))

  return(out)
}

# Primary manuscript analysis:
# full_data <- Read_TSCS2012()
#
# Marriage-stratified analyses:
# married_data <- Read_TSCS2012(marriage_filter = "married")
# nonmarried_data <- Read_TSCS2012(marriage_filter = "nonmarried")
#
# Four-category sensitivity analysis:
# full_data_4cat <- Read_TSCS2012(collapse_z = FALSE)
