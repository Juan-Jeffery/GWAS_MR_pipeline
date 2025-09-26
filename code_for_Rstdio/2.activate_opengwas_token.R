### 嵌入opengwas 令牌 ###
file.edit("~/.Renviron")
#把底下這行放入.Renviron檔案中
#OPENGWAS_JWT=令牌

###restart 令牌
readRenviron("~/.Renviron")

#
library(TwoSampleMR)
library(gwasglue)
library(gwasvcf)
library(VariantAnnotation)
library(MendelianRandomization)
library(dplyr)
library(ggplot2)

# 檢查是否成功 (會看到令牌內容)
Sys.getenv("OPENGWAS_JWT")
ieugwasr::get_opengwas_jwt()
ieugwasr::user()

####################
