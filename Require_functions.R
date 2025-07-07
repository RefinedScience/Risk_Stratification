######################### # Copyright disclaimer and license notice ##################################
######################### # Copyright disclaimer and license notice ##################################
######################### # Copyright disclaimer and license notice ##################################
######################### # Copyright disclaimer and license notice ##################################
######################### # Copyright disclaimer and license notice ##################################
######################### # Copyright disclaimer and license notice ##################################



# ©2025 AML JV, LLC.  All rights reserved.  Use of this software is subject to the AML JV, 
# LLC Non-Commercial Research License. No other use of this software, including, without 
# limitation, use for any commercial purposes, is permitted without the express written consent 
# of AML JV, LLC.



######################### # Copyright disclaimer and license notice ##################################
######################### # Copyright disclaimer and license notice ##################################
######################### # Copyright disclaimer and license notice ##################################
######################### # Copyright disclaimer and license notice ##################################
######################### # Copyright disclaimer and license notice ##################################
######################### # Copyright disclaimer and license notice ##################################



## Definitions of inputs and functions
# udata = Analytical dataset saved as a tibble or dataframe. 


# create new labels for hazard ratio condensing more classes 
numeric_reference_profile_ov <- function(udata){
  
  symbol <- "\u2265"; symbol2 <- intToUtf8(8804)
  
  
  try(udata$DX_LAG_RANGE_V1 <- ifelse(udata$DATE_DIFF <= 5, "[0,5]", 
                                      ifelse(udata$DATE_DIFF > 5, ">5", "Missing")))  
  # date diff between DX and TX start date
  
  # hospitalization
  try(udata$DEM_AGE_RANGE <- ifelse(udata$DEM_AGE <= 75, paste0("", symbol2, "75"), ">75"))
  
  if(sum(is.na(udata$DEM_AGE)) > 0){
    stop("There is a missing age!!!!!")
  }
  
  # hospitalization
  try(udata$HOSP_COUNT_RANGE_V1 <- ifelse(udata$HOSP_COUNT <= 2, "[0,2]", ">2"))
  
  # Transfusions 
  try(udata$LAB_TRANS_RBC_COUNT_RANGE_V1 <- ifelse(udata$LAB_TRANS_RBC_COUNT <= 0, "0", 
                                                   ifelse(udata$LAB_TRANS_RBC_COUNT > 0 & 
                                                            udata$LAB_TRANS_RBC_COUNT <= 5, "[1,5]", ">5")))
  
  try(udata$LAB_TRANS_PLATELET_COUNT_RANGE_V1 <- ifelse(udata$LAB_TRANS_PLATELET_COUNT <= 0, "0", 
                                                        ifelse(udata$LAB_TRANS_PLATELET_COUNT > 0 & 
                                                                 udata$LAB_TRANS_PLATELET_COUNT <= 5, "[1,5]", ">5")))
  
  
  # Blasts
  try(udata$BLASTS_RANGE_V1 <-  ifelse(is.na(udata$BLASTS),  "Missing", 
                                       ifelse(udata$BLASTS <= 25, "[0,25]", ">25")))
  
  
  # new labs thresholds based on lln/uln
  yy <- grep("LAB_", colnames(udata), value = T)
  if(length(grep("_ULN", yy, value = F)) > 1){
    yy <- yy[-grep("_ULN", yy, value = F)]
  }
  if(length(grep("_LLN", yy, value = F)) > 1){
    yy <- yy[-grep("_LLN", yy, value = F)]
  }
  if(length(grep("_TRANS", yy, value = F)) > 1){
    yy <- yy[-grep("_TRANS", yy, value = F)]
  }
  
  if(length(grep("_RANGE_V1", yy, value = F)) >= 1){
    yy <- yy[-grep("_RANGE_V1", yy, value = F)] # needed within imputation as by then - there will be already _RANGE_V1 
  }
  
  if(length(grep("_RANGE_V999", yy, value = F)) >= 1){
    yy <- yy[-grep("_RANGE_V999", yy, value = F)] # needed within imputation as by then - there will be already _RANGE_V1 
  }
  
  if(length(yy) > 0) {
    for(jj in 1 : length(yy)){
      if(yy[jj] %!in% c("LAB_ALT", "LAB_AST", "LAB_BILIRUBIN", "LAB_CREATININE", 
                        "LAB_LDH", "LAB_URICACID")){
        QQ <- paste0("try(udata$", yy[jj],"_RANGE_V1 <- unlist(lapply(seq_len(nrow(udata)), function(ii) {
              ifelse(is.na(udata$", yy[jj],"[ii]), 'Missing',  ifelse(udata$", yy[jj],"[ii] < udata$", 
              yy[jj],"_LLN[ii] | udata$", yy[jj],"[ii] > udata$", yy[jj],"_ULN[ii], '<LLN/>ULN', 'Normal'))
              })))")
        QQ <- eval(parse(text = QQ))
      } else if(yy[jj] %in% c("LAB_ALT", "LAB_AST", "LAB_BILIRUBIN", "LAB_CREATININE", 
                              "LAB_LDH", "LAB_URICACID")){
        QQ <- paste0("try(udata$", yy[jj],"_RANGE_V1 <- unlist(lapply(seq_len(nrow(udata)), function(ii) {
              ifelse(is.na(udata$", yy[jj],"[ii]), 'Missing',  ifelse(udata$", 
              yy[jj],"[ii] > udata$", yy[jj],"_ULN[ii], '>ULN', 'Normal'))
              })))")
        QQ <- eval(parse(text = QQ))
      }    
    }
  }
  
  
  # hgb
  try(udata$LAB_HEMOGLOBIN_RANGE_V1 <-   ifelse(is.na(udata$LAB_HEMOGLOBIN), "Missing",  
                                                ifelse(udata$LAB_HEMOGLOBIN < 8, "<8", paste0("", symbol, "8"))))
  # plat
  try(udata$LAB_PLATELETS_RANGE_V1 <- ifelse(is.na(udata$LAB_PLATELETS), "Missing",  
                                             ifelse(udata$LAB_PLATELETS < 60, "<60", paste0("", symbol, "60") )))
  # wbc
  try(udata$LAB_WBC_RANGE_V1 <-  ifelse(is.na(udata$LAB_WBC), "Missing",  
                                        ifelse(udata$LAB_WBC < 1.1 | udata$LAB_WBC > 11.1, "<1.1|>11.1", "[1.1,11.1]")))
  
  
  # ldh
  try(udata$LAB_LDH_RANGE_V999 <- ifelse(is.na(udata$LAB_LDH), "Missing",  
                                         ifelse(udata$LAB_LDH < 271, "<271", paste0("", symbol, "271"))))
  
  
  # albumin
  try(udata$LAB_ALBUMIN_RANGE_V999 <- ifelse(is.na(udata$LAB_ALBUMIN), "Missing",  
                                             ifelse(udata$LAB_ALBUMIN < 3.5, "<3.5", paste0("", symbol, "3.5"))))
  
  
  # alt (modified)
  try(udata$LAB_ALT_RANGE_V999 <- ifelse(is.na(udata$LAB_ALT), "Missing", 
                                         ifelse(udata$LAB_ALT < 30, "<30", paste0("", symbol, "30"))))
  
  
  # ast (modified)
  try(udata$LAB_AST_RANGE_V999 <-  ifelse(is.na(udata$LAB_AST), "Missing", 
                                          ifelse(udata$LAB_AST < 30, "<30",  paste0("", symbol, "30"))))
  
  
  # anc
  try(udata$LAB_ANC_RANGE_V999 <-  ifelse(is.na(udata$LAB_ANC), "Missing",
                                          ifelse(udata$LAB_ANC >= 0.1 & udata$LAB_ANC <= 10, "[0.1,10]", "<0.1|>10")))
  
  # bilirubin (modified)
  try(udata$LAB_BILIRUBIN_RANGE_V999 <- ifelse(is.na(udata$LAB_BILIRUBIN), "Missing", 
                                               ifelse(udata$LAB_BILIRUBIN < 1, "<1", paste0("", symbol, "1"))))
  
  # calcium 
  try(udata$LAB_CALCIUM_RANGE_V999 <- ifelse(is.na(udata$LAB_CALCIUM), "Missing", 
                                             ifelse(udata$LAB_CALCIUM < 8, "<8", paste0("", symbol, "8"))))
  
  # creatinine (modified)
  try(udata$LAB_CREATININE_RANGE_V999 <-  ifelse(is.na(udata$LAB_CREATININE), "Missing", 
                                                 ifelse(udata$LAB_CREATININE < 1, "<1", paste0("", symbol, "1"))))
  
  
  # fibrinogen
  try(udata$LAB_FIBRINOGEN_RANGE_V999 <- ifelse(is.na(udata$LAB_FIBRINOGEN), "Missing", 
                                                ifelse(udata$LAB_FIBRINOGEN < 150 | udata$LAB_FIBRINOGEN > 400, "<150|>400",  "[150,400]")))
  
  # hgb
  try(udata$LAB_HEMOGLOBIN_RANGE_V999 <-   ifelse(is.na(udata$LAB_HEMOGLOBIN), "Missing",  
                                                  ifelse(udata$LAB_HEMOGLOBIN < 8, "<8", paste0("", symbol, "8"))))
  
  # lymph
  try(udata$LAB_LYMPHOCYTES_RANGE_V999 <- ifelse(is.na(udata$LAB_LYMPHOCYTES), "Missing",  
                                                 ifelse(udata$LAB_LYMPHOCYTES < 0.8 |
                                                          udata$LAB_LYMPHOCYTES > 4.8, "<0.8|>4.8", "[0.8,4.8]")))
  
  # phosphorus (modified)
  try(udata$LAB_PHOSPHORUS_RANGE_V999 <- ifelse(is.na(udata$LAB_PHOSPHORUS), "Missing", 
                                                ifelse(udata$LAB_PHOSPHORUS < 4, "<4", paste0("", symbol, "4"))))
  
  # plat
  try(udata$LAB_PLATELETS_RANGE_V999 <- ifelse(is.na(udata$LAB_PLATELETS), "Missing",  
                                               ifelse(udata$LAB_PLATELETS < 60, "<60", paste0("", symbol, "60") )))
  
  # potassium (modified)
  try(udata$LAB_POTASSIUM_RANGE_V999 <- ifelse(is.na(udata$LAB_POTASSIUM), "Missing",  
                                               ifelse(udata$LAB_POTASSIUM < 3.9, "<3.9", paste0("", symbol, "3.9"))))
  
  
  # uric-acid 
  try(udata$LAB_URICACID_RANGE_V999 <- ifelse(is.na(udata$LAB_URICACID), "Missing",  
                                              ifelse(udata$LAB_URICACID < 7, "<7", paste0("", symbol, "7"))))
  
  # wbc
  try(udata$LAB_WBC_RANGE_V999 <-  ifelse(is.na(udata$LAB_WBC), "Missing",  
                                          ifelse(udata$LAB_WBC < 1.1 | udata$LAB_WBC > 11.1, "<1.1|>11.1", "[1.1,11.1]")))
  
  # ejection
  try(udata$EJECTION_FRACTION_RANGE_V1 <- ifelse(is.na(udata$EJECTION_FRACTION), "Missing",  
                                                 ifelse(udata$EJECTION_FRACTION < 50 | 
                                                          udata$EJECTION_FRACTION > 75, "<50|>75", "[50,75]")))
  
  # BMI
  try(udata$BMI_RANGE_V1 <- ifelse(is.na(udata$BMI), "Missing",  ifelse(udata$BMI < 25, "<25", 
                                                                        paste0("", symbol, "25"))))
  
  
  ## setting references for all numeric 
  udata <- udata %>% mutate_if(is.character, as.factor) 
  try(udata <- udata %>% mutate(HOSP_COUNT_RANGE_V1 = fct_relevel(HOSP_COUNT_RANGE_V1, 
              c("[0,2]", ">2")), HOSP_COUNT_RANGE_V1 = relevel(HOSP_COUNT_RANGE_V1, ref = '[0,2]'))); 
  try(udata <- udata %>% mutate(DX_LAG_RANGE_V1 = fct_relevel(DX_LAG_RANGE_V1, 
              c('[0,5]', '>5', "Missing")), DX_LAG_RANGE_V1 = relevel(DX_LAG_RANGE_V1, ref = '[0,5]')));
  try(udata <- udata %>% mutate(LAB_TRANS_RBC_COUNT_RANGE_V1 = fct_relevel(LAB_TRANS_RBC_COUNT_RANGE_V1, c("0", "[1,5]", ">5")), 
              LAB_TRANS_RBC_COUNT_RANGE_V1 = relevel(LAB_TRANS_RBC_COUNT_RANGE_V1, ref = '0'))); 
  try(udata <- udata %>% mutate(LAB_TRANS_PLATELET_COUNT_RANGE_V1 = fct_relevel(LAB_TRANS_PLATELET_COUNT_RANGE_V1, 
              c("0", "[1,5]", ">5")), LAB_TRANS_PLATELET_COUNT_RANGE_V1 = relevel(LAB_TRANS_PLATELET_COUNT_RANGE_V1, ref = '0')));
  try(udata <- udata %>% mutate(BLASTS_RANGE_V1 = fct_relevel(BLASTS_RANGE_V1, c("[0,25]", ">25", "Missing")), BLASTS_RANGE_V1 = relevel(BLASTS_RANGE_V1, ref = '[0,25]'))); 
  try(udata <- udata %>% mutate(BMI_RANGE_V1 = fct_relevel(BMI_RANGE_V1, c("<25", 
                                                                           paste0("", symbol, "25"), "Missing")), BMI_RANGE_V1 = relevel(BMI_RANGE_V1, ref = '<25')));
  try(udata <- udata %>% mutate(EJECTION_FRACTION_RANGE_V1 = fct_relevel(EJECTION_FRACTION_RANGE_V1, c("[50,75]", "<50|>75", "Missing")), EJECTION_FRACTION_RANGE_V1 = relevel(EJECTION_FRACTION_RANGE_V1, ref = "[50,75]")));
  
  try(udata <- udata %>% mutate(LAB_LDH_RANGE_V999 = fct_relevel(LAB_LDH_RANGE_V999, c("<271", paste0("", symbol, "271"), "Missing")), LAB_LDH_RANGE_V999 = relevel(LAB_LDH_RANGE_V999, ref = '<271'))); 
  try(udata <- udata %>% mutate(LAB_ALBUMIN_RANGE_V999 = fct_relevel(LAB_ALBUMIN_RANGE_V999, c("<3.5", paste0("", symbol, "3.5"), "Missing")), LAB_ALBUMIN_RANGE_V999 = relevel(LAB_ALBUMIN_RANGE_V999, ref = paste0("", symbol, "3.5"))))
  try(udata <- udata %>% mutate(LAB_ALT_RANGE_V999 = fct_relevel(LAB_ALT_RANGE_V999, c("<30", paste0("", symbol, "30"), "Missing")), LAB_ALT_RANGE_V999 = relevel(LAB_ALT_RANGE_V999, ref = '<30'))); 
  try(udata <- udata %>% mutate(LAB_AST_RANGE_V999 = fct_relevel(LAB_AST_RANGE_V999, c("<30", paste0("", symbol, "30"), "Missing")), LAB_AST_RANGE_V999 = relevel(LAB_AST_RANGE_V999, ref = '<30')));
  try(udata <- udata %>% mutate(LAB_ANC_RANGE_V999 = fct_relevel(LAB_ANC_RANGE_V999, c("[0.1,10]", "<0.1|>10", "Missing")), LAB_ANC_RANGE_V999 = relevel(LAB_ANC_RANGE_V999, ref = '[0.1,10]')));
  try(udata <- udata %>% mutate(LAB_BILIRUBIN_RANGE_V999 = fct_relevel(LAB_BILIRUBIN_RANGE_V999, c("<1", paste0("", symbol, "1"), "Missing")), LAB_BILIRUBIN_RANGE_V999 = relevel(LAB_BILIRUBIN_RANGE_V999, ref = '<1'))); 
  try(udata <- udata %>% mutate(LAB_CALCIUM_RANGE_V999 = fct_relevel(LAB_CALCIUM_RANGE_V999, c("<8",  paste0("", symbol, "8"), "Missing")), LAB_CALCIUM_RANGE_V999 = relevel(LAB_CALCIUM_RANGE_V999, ref = '<8')));
  try(udata <- udata %>% mutate(LAB_FIBRINOGEN_RANGE_V999 = fct_relevel(LAB_FIBRINOGEN_RANGE_V999, c("[150,400]", "<150|>400", "Missing")), LAB_FIBRINOGEN_RANGE_V999 = relevel(LAB_FIBRINOGEN_RANGE_V999, ref = '[150,400]'))); 
  try(udata <- udata %>% mutate(LAB_HEMOGLOBIN_RANGE_V999 = fct_relevel(LAB_HEMOGLOBIN_RANGE_V999, c("<8",  paste0("", symbol, "8"), "Missing")), LAB_HEMOGLOBIN_RANGE_V999 = relevel(LAB_HEMOGLOBIN_RANGE_V999, ref =  "<8")));
  try(udata <- udata %>% mutate(LAB_LYMPHOCYTES_RANGE_V999 = fct_relevel(LAB_LYMPHOCYTES_RANGE_V999, c("[0.8,4.8]", "<0.8|>4.8", "Missing")), LAB_LYMPHOCYTES_RANGE_V999 = relevel(LAB_LYMPHOCYTES_RANGE_V999, ref = '[0.8,4.8]')));
  try(udata <- udata %>% mutate(LAB_PHOSPHORUS_RANGE_V999 = fct_relevel(LAB_PHOSPHORUS_RANGE_V999, c("<4", paste0("", symbol, "4"), "Missing")), LAB_PHOSPHORUS_RANGE_V999 = relevel(LAB_PHOSPHORUS_RANGE_V999, ref = '<4')));
  try(udata <- udata %>% mutate(LAB_PLATELETS_RANGE_V999 = fct_relevel(LAB_PLATELETS_RANGE_V999, c("<60", paste0("", symbol, "60"), "Missing")), LAB_PLATELETS_RANGE_V999 = relevel(LAB_PLATELETS_RANGE_V999, ref = paste0("", symbol, "60"))));
  try(udata <- udata %>% mutate(LAB_POTASSIUM_RANGE_V999 = fct_relevel(LAB_POTASSIUM_RANGE_V999, c("<3.9", paste0("", symbol, "3.9"), "Missing")), LAB_POTASSIUM_RANGE_V999 = relevel(LAB_POTASSIUM_RANGE_V999, ref = '<3.9'))); 
  try(udata <- udata %>% mutate(LAB_URICACID_RANGE_V999 = fct_relevel(LAB_URICACID_RANGE_V999, c("<7", paste0("", symbol, "7"), "Missing")), LAB_URICACID_RANGE_V999 = relevel(LAB_URICACID_RANGE_V999, ref = '<7')));
  try(udata <- udata %>% mutate(LAB_CREATININE_RANGE_V999 = fct_relevel(LAB_CREATININE_RANGE_V999, c("<1", paste0("", symbol, "1"), "Missing")), LAB_CREATININE_RANGE_V999 = relevel(LAB_CREATININE_RANGE_V999, ref = '<1')));
  
  try(udata <- udata %>% mutate(LAB_WBC_RANGE_V999 = fct_relevel(LAB_WBC_RANGE_V999, c("[1.1,11.1]", "<1.1|>11.1", "Missing")), LAB_WBC_RANGE_V999 = relevel(LAB_WBC_RANGE_V999, ref = '[1.1,11.1]')));
  
  
  
  
  ## new labs 
  try(udata <- udata %>% mutate(LAB_LDH_RANGE_V1 = fct_relevel(LAB_LDH_RANGE_V1, c('Normal', '>ULN', 'Missing')), LAB_LDH_RANGE_V1 = relevel(LAB_LDH_RANGE_V1, ref = 'Normal'))); 
  try(udata <- udata %>% mutate(LAB_ALBUMIN_RANGE_V1 = fct_relevel(LAB_ALBUMIN_RANGE_V1, c('Normal', '<LLN/>ULN', 'Missing')), LAB_ALBUMIN_RANGE_V1 = relevel(LAB_ALBUMIN_RANGE_V1, ref = 'Normal')));
  try(udata <- udata %>% mutate(LAB_ALT_RANGE_V1 = fct_relevel(LAB_ALT_RANGE_V1, c('Normal', '>ULN', 'Missing')), LAB_ALT_RANGE_V1 = relevel(LAB_ALT_RANGE_V1, ref = 'Normal'))); 
  try(udata <- udata %>% mutate(LAB_AST_RANGE_V1 = fct_relevel(LAB_AST_RANGE_V1, c('Normal', '>ULN', 'Missing')), LAB_AST_RANGE_V1 = relevel(LAB_AST_RANGE_V1, ref = 'Normal')));
  try(udata <- udata %>% mutate(LAB_ANC_RANGE_V1 = fct_relevel(LAB_ANC_RANGE_V1, c('Normal', '<LLN/>ULN', 'Missing')), LAB_ANC_RANGE_V1 = relevel(LAB_ANC_RANGE_V1, ref = 'Normal')));
  try(udata <- udata %>% mutate(LAB_BILIRUBIN_RANGE_V1 = fct_relevel(LAB_BILIRUBIN_RANGE_V1, c('Normal', '>ULN', 'Missing')), LAB_BILIRUBIN_RANGE_V1 = relevel(LAB_BILIRUBIN_RANGE_V1, ref = 'Normal'))); 
  try(udata <- udata %>% mutate(LAB_CALCIUM_RANGE_V1 = fct_relevel(LAB_CALCIUM_RANGE_V1, c('Normal', '<LLN/>ULN', 'Missing')), LAB_CALCIUM_RANGE_V1 = relevel(LAB_CALCIUM_RANGE_V1, ref = 'Normal')));
  try(udata <- udata %>% mutate(LAB_FIBRINOGEN_RANGE_V1 = fct_relevel(LAB_FIBRINOGEN_RANGE_V1, c('Normal', '<LLN/>ULN', 'Missing')), LAB_FIBRINOGEN_RANGE_V1 = relevel(LAB_FIBRINOGEN_RANGE_V1, ref = 'Normal'))); 
  #try(udata <- udata %>% mutate(LAB_HEMOGLOBIN_RANGE_V1 = fct_relevel(LAB_HEMOGLOBIN_RANGE_V1, c('Normal', '<LLN/>ULN', 'Missing')), LAB_HEMOGLOBIN_RANGE_V1 = relevel(LAB_HEMOGLOBIN_RANGE_V1, ref = 'Normal')));
  try(udata <- udata %>% mutate(LAB_HEMOGLOBIN_RANGE_V1 = fct_relevel(LAB_HEMOGLOBIN_RANGE_V1, c("<8",  paste0("", symbol, "8"), "Missing")), LAB_HEMOGLOBIN_RANGE_V1 = relevel(LAB_HEMOGLOBIN_RANGE_V1, ref =  "<8")));
  
  try(udata <- udata %>% mutate(LAB_LYMPHOCYTES_RANGE_V1 = fct_relevel(LAB_LYMPHOCYTES_RANGE_V1, c('Normal', '<LLN/>ULN', 'Missing')), LAB_LYMPHOCYTES_RANGE_V1 = relevel(LAB_LYMPHOCYTES_RANGE_V1, ref = 'Normal')));
  try(udata <- udata %>% mutate(LAB_PHOSPHORUS_RANGE_V1 = fct_relevel(LAB_PHOSPHORUS_RANGE_V1, c('Normal', '<LLN/>ULN', 'Missing')), LAB_PHOSPHORUS_RANGE_V1 = relevel(LAB_PHOSPHORUS_RANGE_V1, ref = 'Normal')));
  
  #try(udata <- udata %>% mutate(LAB_PLATELETS_RANGE_V1 = fct_relevel(LAB_PLATELETS_RANGE_V1, c('Normal', '<LLN/>ULN', 'Missing')), LAB_PLATELETS_RANGE_V1 = relevel(LAB_PLATELETS_RANGE_V1, ref = 'Normal')));
  try(udata <- udata %>% mutate(LAB_PLATELETS_RANGE_V1 = fct_relevel(LAB_PLATELETS_RANGE_V1, c("<60", paste0("", symbol, "60"), "Missing")), LAB_PLATELETS_RANGE_V1 = relevel(LAB_PLATELETS_RANGE_V1, ref = paste0("", symbol, "60"))));
  try(udata <- udata %>% mutate(LAB_POTASSIUM_RANGE_V1 = fct_relevel(LAB_POTASSIUM_RANGE_V1, c('Normal', '<LLN/>ULN', 'Missing')), LAB_POTASSIUM_RANGE_V1 = relevel(LAB_POTASSIUM_RANGE_V1, ref = 'Normal'))); 
  try(udata <- udata %>% mutate(LAB_URICACID_RANGE_V1 = fct_relevel(LAB_URICACID_RANGE_V1, c('Normal', '>ULN', 'Missing')), LAB_URICACID_RANGE_V1 = relevel(LAB_URICACID_RANGE_V1, ref = 'Normal')));
  try(udata <- udata %>% mutate(LAB_CREATININE_RANGE_V1 = fct_relevel(LAB_CREATININE_RANGE_V1, c('Normal', '>ULN', 'Missing')), LAB_CREATININE_RANGE_V1 = relevel(LAB_CREATININE_RANGE_V1, ref = 'Normal')));
  
  #try(udata <- udata %>% mutate(LAB_WBC_RANGE_V1 = fct_relevel(LAB_WBC_RANGE_V1, c('Normal', '<LLN/>ULN', 'Missing')), LAB_WBC_RANGE_V1 = relevel(LAB_WBC_RANGE_V1, ref = 'Normal')));
  try(udata <- udata %>% mutate(LAB_WBC_RANGE_V1 = fct_relevel(LAB_WBC_RANGE_V1, c("[1.1,11.1]", "<1.1|>11.1", "Missing")), LAB_WBC_RANGE_V1 = relevel(LAB_WBC_RANGE_V1, ref = '[1.1,11.1]')));
  
  try(udata <- udata %>% mutate(LAB_HEMOGLOBIN_RANGE_V1 = fct_relevel(LAB_HEMOGLOBIN_RANGE_V1, c('<8', '[8,18]', 'Missing')), LAB_HEMOGLOBIN_RANGE_V1 = relevel(LAB_HEMOGLOBIN_RANGE_V1, ref = '[8,18]')));
  try(udata <- udata %>% mutate(LAB_PLATELETS_RANGE_V1 = fct_relevel(LAB_PLATELETS_RANGE_V1, c('<60', '[60,500]', '>500', 'Missing')), LAB_PLATELETS_RANGE_V1 = relevel(LAB_PLATELETS_RANGE_V1, ref = '[60,500]')));
  
  
  try(udata <- udata %>% mutate(DEM_AGE_RANGE = fct_relevel(DEM_AGE_RANGE, c(paste0("", symbol2, "75"), ">75")), DEM_AGE_RANGE = relevel(DEM_AGE_RANGE, ref = paste0("", symbol2, "75")))); 
  
  output = list(result = udata)
}


# create reference labels for categorical 
categorical_reference_profile_ov <- function(udata){
  
  udata <- as.data.frame(udata)
  
  ## fc 
  if(length(grep("FC_", colnames(udata), value = T)) >= 1){
    QQ <- paste0("", grep("FC_", colnames(udata), value = T), " = fct_relevel(", grep("FC_", colnames(udata), value = T), ", c('No', 'Yes', 'Missing')), ", grep("FC_", colnames(udata), value = T)," = 
         relevel(", grep("FC_", colnames(udata), value = T), ", ref = 'No')", collapse = ", ")
    QQ <- paste0("udata <- udata %>% mutate(", QQ, ")")
    QQ <- eval(parse(text = QQ))
  }
  
  ## IHC 
  if(length(grep("IHC_", colnames(udata), value = T)) >= 1){
    QQ <- paste0("", grep("IHC_", colnames(udata), value = T), " = IHCt_relevel(", grep("IHC_", colnames(udata), value = T), ", c('No', 'Yes', 'Missing')), ", grep("IHC_", colnames(udata), value = T)," = 
         relevel(", grep("IHC_", colnames(udata), value = T), ", ref = 'No')", collapse = ", ")
    QQ <- paste0("udata <- udata %>% mutate(", QQ, ")")
    QQ <- eval(parse(text = QQ))
  }
  
  ## combo variables
  if(length(grep("_COMBO", colnames(udata), value = T)) >= 1){
    
    hg <- grep("_COMBO", colnames(udata), value = T)
    
    hgr <- unlist(lapply(seq_len(length(hg)), function(jj)  grep("ven/aza/cusa", unique(udata[, which(colnames(udata) %in% hg[jj]) ]), value = T)))
    #hgr <- unlist(lapply(seq_len(length(hg)), function(jj)  grep("superc20", unique(udata[, which(colnames(udata) %in% hg[jj]) ]), value = T)))
    
    # Deleting (-ve) if there is any 
    ext <- c(grep("(-ve)", hgr, value = F), grep("Intermediate", hgr, value = F), grep("Favorable", hgr, value = F), grep("LDH ≤ ULN", hgr, value = F))
    if(length(ext) >= 1){
      hgr <- hgr[- ext]
    }
    
    QQ <- paste0("", grep("_COMBO", colnames(udata), value = T)," =  relevel(", grep("_COMBO", colnames(udata), value = T), ", ref = '", hgr, "')", collapse = ", ")
    QQ <- paste0("try(udata <- udata %>% mutate(", QQ, "))")
    QQ <- eval(parse(text = QQ))
  }
  
  ## cyt
  if(length(grep("CYT_", colnames(udata), value = T)) >= 1){
    QQ <- paste0("", grep("CYT_", colnames(udata), value = T), " = fct_relevel(", grep("CYT_", colnames(udata), value = T), ", c('No', 'Yes', 'Missing')), ", grep("CYT_", colnames(udata), value = T)," = 
         relevel(", grep("CYT_", colnames(udata), value = T), ", ref = 'No')", collapse = ", ")
    QQ <- paste0("udata <- udata %>% mutate(", QQ, ")")
    QQ <- eval(parse(text = QQ))
  }
  
  
  ## ngs 
  if(length(grep("NGS_", colnames(udata), value = T)) >= 1){
    QQ <- paste0("", grep("NGS_", colnames(udata), value = T), " = fct_relevel(", grep("NGS_", colnames(udata), value = T), ", c('No', 'Yes', 'Missing')), ", grep("NGS_", colnames(udata), value = T)," = 
         relevel(", grep("NGS_", colnames(udata), value = T), ", ref = 'No')", collapse = ", ")
    QQ <- paste0("udata <- udata %>% mutate(", QQ, ")")
    QQ <- eval(parse(text = QQ))
  }
  
  ## fish 
  if(length(grep("FISH_", colnames(udata), value = T)) >= 1){
    QQ <- paste0("", grep("FISH_", colnames(udata), value = T), " = fct_relevel(", grep("FISH_", colnames(udata), value = T), ", c('No', 'Yes', 'Missing')), ", grep("FISH_", colnames(udata), value = T)," = 
         relevel(", grep("FISH_", colnames(udata), value = T), ", ref = 'No')", collapse = ", ")
    QQ <- paste0("udata <- udata %>% mutate(", QQ, ")")
    QQ <- eval(parse(text = QQ))
  }
  
  ## mut 
  if(length(grep("_MUTATION", colnames(udata), value = T)) >= 1){
    QQ <- paste0("", grep("_MUTATION", colnames(udata), value = T), " = fct_relevel(", grep("_MUTATION", colnames(udata), value = T), ", c('No', 'Yes', 'Missing')), ", grep("_MUTATION", colnames(udata), value = T)," = 
         relevel(", grep("_MUTATION", colnames(udata), value = T), ", ref = 'No')", collapse = ", ")
    QQ <- paste0("udata <- udata %>% mutate(", QQ, ")")
    QQ <- eval(parse(text = QQ))
  }
  
  
  ## hx
  if(length(grep("MEDHX_", colnames(udata), value = T)) >= 1){
    QQ <- paste0("", grep("MEDHX", colnames(udata), value = T), " = fct_relevel(", grep("MEDHX", colnames(udata), value = T), ", c('No', 'Yes', 'Missing')), ", grep("MEDHX", colnames(udata), value = T)," = 
         relevel(", grep("MEDHX", colnames(udata), value = T), ", ref = 'No')", collapse = ", ")
    QQ <- paste0("udata <- udata %>% mutate(", QQ, ")")
    QQ <- eval(parse(text = QQ))
  }
  
  
  ## eln/ecog /gender/race
  try(udata <- udata %>% mutate(DEM_ETHNICITY = fct_relevel(DEM_ETHNICITY, c("Hispanic", "Non-Hispanic", "Missing")), DEM_ETHNICITY = relevel(DEM_ETHNICITY, ref = "Non-Hispanic")), silent = T)
  try(udata <- udata %>% mutate(ECOG = fct_relevel(ECOG, c("0", "1", paste0("", symbol, "2"), "Missing")), ECOG = relevel(ECOG, ref = "0")), silent = T)
  
  try(udata <- udata %>% mutate(MOD_RISK = relevel(MOD_RISK, ref = "Favorable")), silent = T)
  try(udata <- udata %>% mutate(TX_TYPE = relevel(TX_TYPE, ref = "ven/aza")), silent = T)
  #try(udata <- udata %>% mutate(TX_TYPE = relevel(TX_TYPE, ref = "ven/aza/cusa")), silent = T) # for VAC vs VA paper 
  #try(udata <- udata %>% mutate(TX_TYPE = relevel(TX_TYPE, ref = "superc20")), silent = T) # for superc20 ref
  
  try(udata <- udata %>% mutate(TX_TYPE2 = relevel(TX_TYPE2, ref = "ven/aza")), silent = T) #/cusa
  try(udata <- udata %>% mutate(ELN_RISK_GROUP = fct_relevel(ELN_RISK_GROUP, c("Favorable", "Intermediate", "Adverse")), ELN_RISK_GROUP = relevel(ELN_RISK_GROUP, ref = "Favorable")), silent = T)
  try(udata <- udata %>% mutate(ELN24_REFINED_RISK_GROUP = fct_relevel(ELN24_REFINED_RISK_GROUP, c("Favorable", "Intermediate", "Adverse")), ELN_RISK_GROUP = relevel(ELN_RISK_GROUP, ref = "Favorable")), silent = T)
  try(udata <- udata %>% mutate(DEM_SEX = fct_relevel(DEM_SEX, c("Male", "Female")), DEM_SEX = relevel(DEM_SEX, ref = "Male")), silent = T)
  try(udata <- udata %>% mutate(DEM_RACE = fct_relevel(DEM_RACE, c("Black", "White", "Missing")), DEM_RACE = relevel(DEM_RACE, ref = "White")), silent = T)
  try(udata <- udata %>% mutate(DEM_AGE_RANGE = fct_relevel(DEM_AGE_RANGE, c(paste0("", symbol2, "65"), "(65,75]", ">75")), DEM_AGE_RANGE = relevel(DEM_AGE_RANGE, ref = paste0("", symbol2, "65"))), silent = T) 
  try(udata <- udata %>% mutate(HOSP_LOS_RANGE_V1 = fct_relevel(HOSP_LOS_RANGE_V1, c('[0,5]', '>5')), HOSP_LOS_RANGE_V1 = relevel(HOSP_LOS_RANGE_V1, ref = '[0,5]')), silent = T);
  try(udata <- udata %>% mutate(HOSP_LOSICU_RANGE_V1 = fct_relevel(HOSP_LOSICU_RANGE_V1, c('[0,5]/non-ICU', '>5/ICU')), HOSP_LOSICU_RANGE_V1 = relevel(HOSP_LOSICU_RANGE_V1, ref = '[0,5]/non-ICU')), silent = T);
  
  # random1
  try(udata <- udata %>% mutate(RANDOM_VAR1 = fct_relevel(RANDOM_VAR1, c('Yes', 'No')), RANDOM_VAR1 = relevel(RANDOM_VAR1, ref = 'No')), silent = T);
  # random1
  try(udata <- udata %>% mutate(RANDOM_VAR2 = fct_relevel(RANDOM_VAR2, c('Yes', 'No')), RANDOM_VAR2 = relevel(RANDOM_VAR2, ref = 'No')), silent = T);
  # random1
  try(udata <- udata %>% mutate(RANDOM_VAR3 = fct_relevel(RANDOM_VAR3, c('Yes', 'No')), RANDOM_VAR3 = relevel(RANDOM_VAR3, ref = 'No')), silent = T);
  # random1
  try(udata <- udata %>% mutate(RANDOM_VAR4 = fct_relevel(RANDOM_VAR4, c('Yes', 'No')), RANDOM_VAR4 = relevel(RANDOM_VAR4, ref = 'No')), silent = T);
  # random1
  try(udata <- udata %>% mutate(RANDOM_VAR5 = fct_relevel(RANDOM_VAR5, c('Yes', 'No')), RANDOM_VAR5 = relevel(RANDOM_VAR5, ref = 'No')), silent = T);
  
  
  try(udata <- udata %>% mutate(CD7_PRIM_CATG = fct_relevel(CD7_PRIM_CATG, c(paste0("", symbol2, "5"), "(5,15]", "(15,25]", ">25")), 
                                CD7_PRIM_CATG = relevel(CD7_PRIM_CATG, ref = c(paste0("", symbol2, "5")))), silent = T);
  
  try(udata <- udata %>% mutate(CD7_MONO_CATG = fct_relevel(CD7_MONO_CATG, c(paste0("", symbol2, "5"), "(5,15]", "(15,25]", ">25")), 
                                CD7_MONO_CATG = relevel(CD7_MONO_CATG, ref = c(paste0("", symbol2, "5")))), silent = T);
  
  try(udata <- udata %>% mutate(CD7_BLASTS_CATG = fct_relevel(CD7_BLASTS_CATG, c(paste0("", symbol2, "5"), "(5,15]", "(15,25]", ">25")), 
                                CD7_BLASTS_CATG = relevel(CD7_BLASTS_CATG, ref = c(paste0("", symbol2, "5")))), silent = T);
  
  try(udata <- udata %>% mutate(CD70_BLASTS_CATG = fct_relevel(CD70_BLASTS_CATG, c(paste0("", symbol2, "10"), ">10", "Missing")), 
                                CD70_BLASTS_CATG = relevel(CD70_BLASTS_CATG, ref = paste0("", symbol2, "10")     )), silent = T);
  
  output = list(result = udata)
}


## Definitions of inputs 
# udata = Analytical dataset saved as a tibble or dataframe. 
# impute = Select methods for inverse probability weight or complete cases
# var_to_add = features to add in adjustment models for missing data methods
# only_missing = a Boolean operator taking values TRUE or FALSE (default); FALSE corresponds to full analytically set treating missing element as category 
whatnow_resp_paperML <- function(udata, impute = "MICE", var_to_add, only_missing = FALSE) {
  
  if(only_missing == FALSE){
    
    Main <- udata
    
    if(length(which(Main$time == 0))){ # adjustment if last known alive date is as same as CR/CRi/etc date leading to 0 followup/censoring time for survival
      Main$time[which(Main$time == 0)] <- Main$time[which(Main$time == 0)] + 1 
    }
    
    ## Generate KM for ven/aza data 
    symbol <- "\u2265"; symbol2 <- intToUtf8(8804)
    
    Main <- Main %>% mutate(DEM_AGE_RANGE = ifelse(Main$DEM_AGE <= 75, paste0("", symbol2, "75"), 
                                                   ifelse(Main$DEM_AGE > 75, ">75", NA)))
    Main <- Main %>% mutate_if(is.character, as.factor) %>% mutate(DEM_AGE_RANGE = relevel(DEM_AGE_RANGE, ref = paste0("", symbol2, "75")))
    
    if("Yes" %in% Main$SCT_FLAG | "Y" %in% Main$SCT_FLAG){
      KMdata <- Main %>% mutate(
        VA_FLAG = ifelse(Main$SCT_FLAG == "Yes" & Main$SEC73_FLAG == "No", "ven/aza + SCT", 
                         ifelse(Main$SCT_FLAG == "Yes" & Main$SEC73_FLAG == "Yes", "ven/aza + 7/3 + SCT", 
                                ifelse(Main$SCT_FLAG == "No" & Main$SEC73_FLAG == "Yes", "ven/aza + 7/3", "ven/aza" ))), 
        SCT_FLAG = ifelse(Main$SCT_FLAG == "Yes", "ven/aza + SCT", "ven/aza"), 
        SEC73_FLAG = ifelse(Main$SEC73_FLAG == "Yes", "ven/aza + 7/3", "ven/aza"))
      
      KMdata <- KMdata %>% mutate_if(is.character, as.factor) %>% mutate(SCT_FLAG = relevel(SCT_FLAG, ref = "ven/aza"), 
                                                                         SEC73_FLAG = relevel(SEC73_FLAG, ref = "ven/aza"), 
                                                                         VA_FLAG = relevel(VA_FLAG, ref = "ven/aza")
                                                                      
      )
    } else{
      KMdata <- Main
    }
    sum(is.na(KMdata$time))
  } else if(only_missing == TRUE){
    KMdata <- Main
  }
  
  if(length(which(KMdata$time == 0)) >= 1){
    KMdata$time[which(KMdata$time == 0)] <- 1 # if 0 add 1 on time
  }
  mis_mrn <- NULL
  
  KMdata$missing_wt <- rep(1, nrow(KMdata))
  
  ## Method 4: MICE
  if(impute == "MICE"){
    nm <- unlist(lapply(seq_len(length(var_to_add)), function(jj) {
      QQ <- paste0("length(which(KMdata$", var_to_add[jj], " == 'Missing'))")
      QQ <- eval(parse(text = QQ))
      QQ
    })) 
    nm
    if(length (which(nm >= length(which(is.na(KMdata$censor) == FALSE)) / 3 )) > 0) {
      var_to_add1 <- var_to_add[- which(nm >= length(which(is.na(KMdata$censor) == FALSE)) / 3 )]
      var_to_add <- var_to_add1
    }
    length(var_to_add)
  } else if(impute == "DELETE_CASE"){
    
    ## Method 2: listwise deletion 
    ## Number of missing values for each variable 
    nm <- unlist(lapply(seq_len(length(var_to_add)), function(jj) {
      QQ <- paste0("length(which(KMdata$", var_to_add[jj], " == 'Missing'))")
      QQ <- eval(parse(text = QQ))
      QQ
    })) 
    nm
    if(length(which(nm >= round(nrow(KMdata) / 10))) >= 1){
      var_to_add1 <- var_to_add[-which(nm >= round(nrow(KMdata) / 10))]
    } else if(length(which(nm >= round(nrow(KMdata) / 10))) == 0){
      var_to_add1 <- var_to_add[-which(nm >= 30)]
    }
    
    length(var_to_add1)
    
    ## Find which MRNs are missing 
    mis_mrn <- unique(unlist(lapply(seq_len(length(var_to_add1)), function(jj) {
      QQ <- paste0("KMdata$MRN[which(KMdata$", var_to_add1[jj], " == 'Missing')]")
      QQ <- eval(parse(text = QQ))
      QQ
    }))) 
    length(mis_mrn)
    
    
    
    # the following only for MRD 
    udata = KMdata %>% filter(!is.na(censor)) %>% subset(MRN %!in% mis_mrn)
    unq <- unlist(lapply(seq_len(length(var_to_add1)), function(jj) {
      QQ <- paste0("length(unique(udata$", var_to_add1[jj], "))")
      QQ <- eval(parse(text = QQ))
      QQ
    })) 
    which(unq  == 1)
    if(length(which(unq  == 1)) >= 1){
      var_to_add1 <- var_to_add1[-which(unq  == 1)]
    }
    var_to_add <- var_to_add1
    
  } else if(impute == "IPW" & resp_type %in% c("FUP30D_RESP", "OPTIMUM_RESP_COMB")){
    
    
    ## Method 2: listwise deletion 
    ## Number of missing values for each variable 
    nm <- unlist(lapply(seq_len(length(var_to_add)), function(jj) {
      QQ <- paste0("length(which(KMdata$", var_to_add[jj], " == 'Missing'))")
      QQ <- eval(parse(text = QQ))
      QQ
    })) 
    nm
    
    # first we dlete pts with 1/7th of total data is missing (~15%). If we don't find any variable like that then we delete variables with at-least 30 missing cases 
    if(length(which(nm >= round(nrow(KMdata) / 10))) >= 1){ # if volume of missing data on a var is greater than 1/7th (~15%) of total N
      var_to_add1 <- var_to_add[-which(nm >= round(nrow(KMdata) / 10))]
    } else if(length(which(nm >= round(nrow(KMdata) / 10))) == 0){ 
      if(length(which(nm >= 30)) >= 1){  #if volume of missing data on a var is NOT greater than 1/7th of total N BUT greater than >30
        var_to_add1 <- var_to_add[-which(nm >= 30)] # [We can basically say ]
      } else{
        var_to_add1 <- var_to_add
      }
    }
    length(var_to_add1)
    
    if(length(var_to_add1) == 1){
      var_to_add1 <- c(var_to_add1, "DEM_SEX") # just adding a variable to run the model (this potentially a dummy model without any variables)
    }
    
    ## Find which MRNs are missing 
    mis_mrn <- unique(unlist(lapply(seq_len(length(var_to_add1)), function(jj) {
      QQ <- paste0("KMdata$MRN[which(KMdata$", var_to_add1[jj], " == 'Missing')]")
      QQ <- eval(parse(text = QQ))
      QQ
    }))) 
    length(mis_mrn)
    
    
    # the following only for MRD 
    udata = KMdata %>% filter(!is.na(censor)) %>% subset(MRN %!in% mis_mrn)
    unq <- unlist(lapply(seq_len(length(var_to_add1)), function(jj) {
      QQ <- paste0("length(unique(udata$", var_to_add1[jj], "))")
      QQ <- eval(parse(text = QQ))
      QQ
    })) 
    which(unq  == 1)
    if(length(which(unq  == 1)) >= 1){
      var_to_add1 <- var_to_add1[-which(unq  == 1)]
    }
    var_to_add <- var_to_add1
    
    
    ## Method 3: IPW (for FUP30D)
    KMdata$msg_status <- rep(1, nrow(KMdata)) 
    KMdata$msg_status[which(KMdata$MRN %in% mis_mrn)] <- 0
    
    if(resp_type == "FUP30D_RESP"){
      rr <- KMdata$FUP30D_RESP
      r1 <- rr[which(is.na(rr) == FALSE)]
    } else if(resp_type == "OPTIMUM_RESP_COMB"){
      rr <- KMdata$OPTIMUM_RESP_COMB
      r1 <- rr[which(is.na(rr) == FALSE)]
    } 
    
    ## adding variables into model along with response
    sbv <- KMdata[which(is.na(rr) == FALSE), ]
    sbv$r1 <- r1
    
    ## checking if it has at-least 4 cases in each medhx variable 
    sbv1 <- sbv %>% select(starts_with("MEDHX_"))
    delhx <- colnames(sbv1)[which(unlist(lapply(seq_len(ncol(sbv1)), function(ii) {
      ifelse(length(which(sbv1[, ii] == "Yes")) >= 4, 1, 0)
    })) == 0)]
    
    if(length(delhx) > 0){
      sbv <- sbv[, colnames(sbv) %!in% delhx] 
    }
    
    ot <- c("DEM_SEX", "DEM_AGE", grep("MEDHX_", colnames(sbv), value = T))
    a1 <- paste0("", ot, "", collapse = " + ")
    a2 <- paste0("r1 + ", a1, " ")
    
    QQ <- paste0("fit_wtg <- glm(as.factor(msg_status) ~ ", a2,", data = sbv,  family = 'binomial')")
    QQ <- eval(parse(text = QQ))
    
    KMdata$missing_wt[which(is.na(KMdata$FUP30D_RESP) == FALSE)] <- 1 / fit_wtg$fitted.values
    
  } else if(impute == "IPW" & resp_type %!in% c("FUP30D_RESP", "OPTIMUM_RESP_COMB")){
    
    ## Method 2: listwise deletion 
    ## Number of missing values for each variable 
    nm <- unlist(lapply(seq_len(length(var_to_add)), function(jj) {
      QQ <- paste0("length(which(KMdata$", var_to_add[jj], " == 'Missing'))")
      QQ <- eval(parse(text = QQ))
      QQ
    })) 
    nm
    
    # first we dlete pts with 1/7th of total data is missing (~15%). If we don't find any variable like that then we delete variables with at-least 30 missing cases 
    if(length(which(nm >= round(nrow(KMdata) / 10))) >= 1){ # if volume of missing data on a var is greater than 1/10th (~10%) of total N
      var_to_add1 <- var_to_add[-which(nm >= round(nrow(KMdata) / 10))]
    } else if(length(which(nm >= round(nrow(KMdata) / 10))) == 0){ 
      if(length(which(nm >= 30)) >= 1){  #if volume of missing data on a var is NOT greater than 1/10th of total N BUT greater than >30
        var_to_add1 <- var_to_add[-which(nm >= 30)] # [We can basically say ]
      } else{
        var_to_add1 <- var_to_add
      }
    }
    
    if(length(var_to_add1) == 1){
      
      var_to_add1 <- c(var_to_add1, "DEM_SEX") # just adding a variable to run the model (this potentially a dummy model without any variables)
    }
    
    ## Find which MRNs are missing 
    mis_mrn <- unique(unlist(lapply(seq_len(length(var_to_add1)), function(jj) {
      QQ <- paste0("KMdata$MRN[which(KMdata$", var_to_add1[jj], " == 'Missing')]")
      QQ <- eval(parse(text = QQ))
      QQ
    }))) 
    length(mis_mrn)
    
    
    # the following only for MRD 
    udata = KMdata %>% filter(!is.na(censor)) %>% subset(MRN %!in% mis_mrn)
    unq <- unlist(lapply(seq_len(length(var_to_add1)), function(jj) {
      QQ <- paste0("length(unique(udata$", var_to_add1[jj], "))")
      QQ <- eval(parse(text = QQ))
      QQ
    })) 
    which(unq  == 1)
    if(length(which(unq  == 1)) >= 1){
      var_to_add1 <- var_to_add1[-which(unq  == 1)]
    }
    var_to_add <- var_to_add1
    
    
    # Method 3: IPW (for censor)
    KMdata$msg_status <- rep(1, nrow(KMdata)) 
    KMdata$msg_status[which(KMdata$MRN %in% mis_mrn)] <- 0 # weight wi used is the inverse of the probability that individual i is a complete case.
    table(KMdata$msg_status)
    
    ## adding variables into model along with response
    sbv <- KMdata %>% filter(!is.na(censor))
    
    ## checking if it has at-least 4 cases in each medhx variable 
    sbv1 <- sbv %>% select(starts_with("MEDHX_"))
    delhx <- colnames(sbv1)[which(unlist(lapply(seq_len(ncol(sbv1)), function(ii) {
      ifelse(length(which(sbv1[, ii] == "Yes")) >= 4, 1, 0)
    })) == 0)]
    
    if(length(delhx) > 0){
      sbv <- sbv[, colnames(sbv) %!in% delhx] 
    }
    
    ot <- c("DEM_SEX", "DEM_AGE", grep("MEDHX_", colnames(sbv), value = T))
    a1 <- paste0("", ot, "", collapse = " + ")
    a2 <- paste0("censor + time + ", a1, " ")
    
    QQ <- paste0("fit_wtg <- glm(as.factor(msg_status) ~ ", a2,", data = sbv,  family = 'binomial')")
    QQ <- eval(parse(text = QQ))
    
    #fit_wtg <- glm(as.factor(msg_status) ~ censor + time, data = , family = "binomial")
    KMdata$missing_wt[which(is.na(KMdata$censor) == FALSE)] <- 1 / fit_wtg$fitted.values
    summary(KMdata$missing_wt)
    
    
  }
  
  output = list(KMdata = KMdata, takem0 = var_to_add, mis_mrn = mis_mrn)
}




## fitting multivariate models
# udata =   Analytical dataset saved as a tibble or dataframe.
# tx_type = Indicator for treatment and control group if we need to evaluate by different treatment combinations (not needed for non-causal model)
# method =  "CoxPH + Penalty + Boot" or "GLM + Penalty + Boot"
# var_to_add = Features to be added in multivariate model
# what = Label set for reference 
# alpha = Significance level
# delete_np = Boolean operator with TRUE or FALSE whether missing data needs to be deleted
# delete_NP_param = Boolean operator whether parameters associated with missing category needs to be excluded; default is TRUE as primary interest is to assess the relationship between Yes/No and response
# logistic = Boolean operator whether the logistic or hazard model needs to be fitted
# nnt = Boolean operator whether marginal risk needs to be estimated
# psm = Boolean operator whether propensity score models need to be estimated (not needed for non-causal model)
# psm_method = Type of propensity score methods (not needed for non-causal model)
# psm_model = Type of models to estimate propensity scores (not needed for non-causal model)
# psm_var = Features to adjust in propensity score model (not needed for non-causal model)
# alltx = Type of treatments needed to estimate for estimating marginal risk (not needed for non-causal model)
# response = Type of follow-up response (30 refers to 30-day response)
# resp_type = "OPTIMUM_RESP_WH_COMB_14d" for best-response and "OS" for overall survival
# B = Number of bootstrap runs
# ngamma = Number of values in the grid for regularization parameter for elastic net
# relax_me = Boolean operator indicating whether elastic net needs to be fitted. FALSE corresponds to L2 norm (ridge) penalty. 
# AA = Seed number
# nfolds = Number of cross-validation
# minnum = Minimum number of cases set for each cross-validation (if it is not null, then we select dynamically the number of cross-validation)
# alpha_val = Values to indicate whether to run LASSO/Ridge/Elastic net penalty. Default is 0 indicating L2 norm penalty.  
# measure_type = Minimization/Maximization metrics to evaluate model performance 
# family_type = "Cox" for survival  or "Binomial" for binary outcome
# imputeme = Boolean operator whether MICE needs to run to impute values for var_to_add features 
# nsample = Number of bootstrap samples 
Inter_TX_multivariate_profile_paperML <- function(udata, tx_type, method, var_to_add,
                                             what = "N", 
                                             alpha = 0.05,
                                             delete_np = FALSE,
                                             delete_NP_param = TRUE,
                                             logistic = FALSE,
                                            
                                             nnt = FALSE, 
                                             
                                             psm = FALSE, 
                                             psm_method = "IPTW", 
                                             
                                             psm_model= NULL, 
                                             psm_var = c("DEM_AGE", "ECOG"),
                                             alltx = "drop vc", 
                                             response = 30,
                                             resp_type = "OS",
                                             B = 100,
                                             nfolds = 15, 
                                             ngamma = 30, 
                                             relax_me = FALSE, 
                                             AA = 1200,
                                             minnum = 5,
                                             alpha_val = 0,
                                             measure_type = "deviance",
                                             family_type = "cox",
                                             imputeme = FALSE,
                                             nsample = 5) {
  
  udata <- as.data.frame(udata)
  
  ## The following will help fitting one feature set for imputed model and another for outcome model
  if(is.list(var_to_add) == TRUE){
    var_to_impute <- var_to_add[[2]]
    var_to_add <- var_to_add[[1]]
  } else{
    var_to_impute <- var_to_add
  }
  
  
  symbol <- "\u2265"; symbol2 <- intToUtf8(8804)
  ME = abs(qnorm(alpha / 2))
  
  ## deleting NPs and its associated subject for better interpretation
  udata <- udata %>% mutate_if(is.factor, as.character) # deleting original factoring labels
  
  if(delete_np == TRUE){
    udata <- udata %>% na_if("NP") %>% na_if("Not performed") %>% na_if("Unable to assess") %>% na_if("Unknown")  %>% na_if("Missing")
    ## deleting patients for NPs
    adj2 <- which(as.vector(apply(udata, 1, function(ff) sum(is.na(ff)))) >= 1) # row adjustment (subjects)
    if(length(adj2) >= 1){
      udata <- udata[- adj2, ]
    }
  }
  
  
  ## Categorizing numerical lab variables
  if(length(grep("RANGE", var_to_add, value = T)) > 0){ 
   try(cval <- numeric_reference_profile_ov(udata)$result)
  }
  
  
  finq025 <- function(ff) {quantile(ff, 0.025)}
  finq975 <- function(ff) {quantile(ff, 0.975)}
  
  ## checking if each catgeory (except missing) label in var_to_add has at-least minn values 
  notry <- grep("_COMBO", var_to_add, value = T)
  dll <- which(
    unlist(lapply(seq_len(length(var_to_add)), function(jj) {
      
      QQ <- paste0("table(cval$", var_to_add[jj], ")")
      QQ <- eval(parse(text = QQ))
      
      if(var_to_add[jj] %!in% notry & length(which(QQ[names(QQ) %!in% "Missing"] < minnum)) >= 1){
        1111 # delete me 
      } else if(var_to_add[jj] %!in% notry & length(QQ) == 1) { # just one label (all Y/ all N/ all missing)
        1111 # delete me 
      } else{
        0
      }
      
    })) == 1111)
  if(length(dll) >= 1){
    var_to_add <- var_to_add[- dll] # delete these and update var_to_add 
    var_to_impute <- var_to_impute[- dll]
  }
  
  
  cval <- cval[, colnames(cval) %in% c("MRN", "time", "censor", 
                                       "TX_TYPE", 
                                       "FUP30D_RESP", 
                                       "OPTIMUM_RESP_WOH_COMB", "OPTIMUM_RESP_WH_COMB",
                                       "OPTIMUM_RESP_WOH_COMB_14d", "OPTIMUM_RESP_WH_COMB_14d", 
                                       var_to_add, var_to_impute, resp_type, psm_var, 
                                       "missing_wt", "ate_wt1")]
  
  dim(cval)
  
  ## convert factor into character
  cval <- cval %>% mutate_if(is.character, as.factor) # re-creating factors
  try(cval <-  categorical_reference_profile_ov(udata = cval)$result)
  
  
  ## If imputed models are needed & whether propensity score models are needed to be fitted 
  if(imputeme == TRUE & psm == FALSE){  
    system.time(mc <- multivariate_impute_ov(udata, varb = var_to_impute, delete_np, resp_type, nsample = nsample, AA = AA, psm = FALSE, psm_var = psm_var, alltx = alltx, psm_model= psm_model))
  } else if(imputeme == TRUE & psm == TRUE & is.null(psm_model) == FALSE){ # Generated psm weights within MICE step for each imputed datatset 
    system.time(mc <- multivariate_impute_ov(udata, varb = var_to_impute, delete_np, resp_type, nsample = nsample, AA = AA, psm = TRUE, psm_var = psm_var, alltx = alltx, psm_model= psm_model))
  }
  
  
  ## Logistic regression
  vifm <- dist_select <- ph_check <- NULL
  
  if(length(var_to_add) == 1){
    if(var_to_add == "TX_TYPE"){
      var_to_add <- c("VAR1", "TX_TYPE")
      cval$VAR1 <- rep(1, nrow(cval))
    }
  }
  
  ## logistic witout imputation
  if(logistic == TRUE & imputeme == FALSE){
    
    nasum <- function(xx) {sum(is.na(xx))}
   

    if(resp_type == "FUP30D_RESP"){ # for 30d response
      cval <- cval %>% filter(!is.na(FUP30D_RESP))
      resp <- cval$FUP30D_RESP
    } else if(resp_type == "OPTIMUM_RESP_WOH_COMB"){ # Best response without CRh 
      cval <- cval %>% filter(!is.na(OPTIMUM_RESP_WOH_COMB))
      resp <- cval$OPTIMUM_RESP_WOH_COMB
    } else if(resp_type == "OPTIMUM_RESP_WH_COMB"){ # Best response with CRh
      cval <- cval %>% filter(!is.na(OPTIMUM_RESP_WH_COMB))
      resp <- cval$OPTIMUM_RESP_WH_COMB
    } else if(resp_type == "OPTIMUM_RESP_WOH_COMB_14d"){ # Best response with without CRh having novel definition
      cval <- cval %>% filter(!is.na(OPTIMUM_RESP_WOH_COMB_14d))
      resp <- cval$OPTIMUM_RESP_WOH_COMB_14d
    } else if(resp_type == "OPTIMUM_RESP_WH_COMB_14d"){ # Best response with CRh having novel definition
      cval <- cval %>% filter(!is.na(OPTIMUM_RESP_WH_COMB_14d))
      resp <- cval$OPTIMUM_RESP_WH_COMB_14d
    } else{
      stop("No response selected")
    }
    
    a1 <- paste0("glm(resp ~ ")
    a2 <- paste0("", var_to_add, "", collapse = " + ")
    
    ww <- paste0("", a1, "", a2, ", data = cval, family = binomial(link = 'logit'), maxit = 200)") # remove NAs by default
    fom <- try(eval(parse(text = ww)))
    
    if(method == "GLM + Penalty + Boot"){
      
      # extract model matrix
      mf <- model.matrix(fom)[, -1]
      my <- resp
      
      
      set.seed(AA + 9999) 
      
      # clipping weights 
      cval$missing_wt <- ifelse(cval$missing_wt >  quantile(cval$missing_wt, 0.95, na.rm = T),  quantile(cval$missing_wt, 0.95, na.rm = T), cval$missing_wt) # clipped at 95% weight value
      
      
      system.time(cv.fit <- cv.glmnet(mf, my, weights = cval$missing_wt, family = family_type, nfolds = ifelse(is.null(minnum) == FALSE, ceiling(nrow(cval) / minnum), 
                                                                                                               nfolds), type.measure = measure_type, #"class", #"auc", #"class", 
                                      gamma = seq(0, 1, length.out = ngamma), alpha = alpha_val, relax = relax_me, parallel = FALSE))
      beta.fit <- coef(cv.fit, s = "lambda.min") 
      name_coef <- rownames(beta.fit)
     
      cores <- detectCores()
      cl <- makeCluster(cores[1] - 2) #not to overload your computer
      registerDoParallel(cl)
      
      system.time(
        betaB <- foreach(bb = 1 : B, .combine = cbind, .packages = c('glmnet',  'ranger', 'doParallel', 'survival')) %dopar% {
         
          set.seed(bb + 8888 + AA)
          n <- nrow(cval)
          wt <- rexp(n, 1)
          wt <- wt / sum(wt)
         
          system.time(cv.fitB <- cv.glmnet(mf, my, family = family_type, weights = cval$missing_wt * wt, 
                                           nfolds = ifelse(is.null(minnum) == FALSE, ceiling(nrow(cval) / minnum), nfolds), 
                                           type.measure = measure_type, alpha = alpha_val,   
                                           gamma = seq(0, 1, length.out = ngamma), relax = relax_me, parallel = FALSE))
          coef(cv.fitB, s = "lambda.min") 
        })
      
      
      stopCluster(cl)
      rn <- rownames(beta.fit)
      
      ## bias-adjusted on fractional-random weight bootstrap https://par.nsf.gov/servlets/purl/10155761
      boot_bound <- do.call(rbind, lapply(seq_len(length(beta.fit)), function(jj) {
      
        ic <- which(rn == rn[jj]) 
        icb <- which(rownames(betaB) == rn[jj])
        bb <- sort(betaB[icb, ])
        q <- length(which(bb < beta.fit[ic])) / B
        lb_ix <- round(B * pnorm((2 * qnorm(q, 0, 1)) + qnorm(alpha / 2, 0, 1)))
        lb_ix <- ifelse(lb_ix == 0, 1, lb_ix)
        ub_ix <- round(B * pnorm((2 * qnorm(q, 0, 1)) + qnorm(1 - (alpha / 2), 0, 1)))
        
        data.frame("names" =  rn[jj], "lb_bias2" = bb[lb_ix], "ub_bias2" = ifelse(ub_ix != 0, bb[ub_ix], quantile(bb, 0.975)))
        
      }))
      
      bias <- (2 * as.numeric(beta.fit)) - apply(betaB, 1, mean) 
      dat <- data.frame("Variable" = name_coef,
                        "bOR" = exp(as.numeric(beta.fit)),
                        "lb_025" = exp(apply(betaB, 1, finq025)),
                        "ub_975" = exp(apply(betaB, 1, finq975)),
                        "sd" = sqrt(apply(betaB, 1, var)),
                        
                        "p_val" =  2 * (1 - pnorm(abs(as.numeric(beta.fit) / sqrt(apply(betaB, 1, var))))),
                        
                        "lb_bias" = exp(bias - ME * sqrt(apply(betaB, 1, var)) ),
                        "ub_bias" = exp(bias + ME * sqrt(apply(betaB, 1, var)) ),
                        
                        "lb_bias2" = exp(boot_bound$lb_bias2), 
                        "ub_bias2" = exp(boot_bound$ub_bias2),
                        
                        "lb_emp" = exp(2 * as.numeric(beta.fit) - apply(betaB, 1, finq975)),
                        "ub_emp" = exp(2 * as.numeric(beta.fit) + apply(betaB, 1, finq025)),
                       
                        "lb" = exp(as.numeric(beta.fit) - ME * sqrt(apply(betaB, 1, var)) ),  
                        "ub" = exp(as.numeric(beta.fit) + ME * sqrt(apply(betaB, 1, var)) ))
      
    }
    
  }
  
  # logistic with imputation
  if(logistic == TRUE & imputeme == TRUE){
    
    dat_imp <- bm <- beta.fit <- list(); 
    
    for(ll in 1 : nsample){ # run for each imputed dataset 
      
      var_to_add_adj <- var_to_add[var_to_add %!in% mc$dropme]
      cval_imp <- mc$subd_imp[[ll]][, colnames(mc$subd_imp[[ll]]) %!in% mc$dropme]
      
      
      nasum <- function(xx) {sum(is.na(xx))}
      apply(cval_imp, 2, nasum)
      
      
      if(length(which(apply(cval_imp[, colnames(cval_imp) %!in% c(grep("FUP30D_", colnames(cval_imp), value = T),
                                                                  grep("FC_CD70", colnames(cval_imp), value = T),
                                                                  grep("CD7_PRIM", colnames(cval), value = T),
                                                                  grep("CD7_MONO", colnames(cval), value = T),
                                                                  grep("CD7_BLASTS", colnames(cval), value = T),
                                                                  grep("ate_wt", colnames(cval_imp), value = T),
                                                                  grep("_COMBO", colnames(cval), value = T),
                                                                  grep("l_", colnames(cval_imp), value = T),
                                                                  grep("OPTIMUM_", colnames(cval_imp), value = T))], 2, nasum) >= 1)) >= 1){
        stop("Missing labs. Add unknown label")
      }
      
      
      if(resp_type == "FUP30D_RESP"){
        cval_imp <- cval_imp %>% filter(!is.na(FUP30D_RESP))
        resp <- cval_imp$FUP30D_RESP
      } else if(resp_type == "OPTIMUM_RESP_WOH_COMB"){
        cval_imp <- cval_imp %>% filter(!is.na(OPTIMUM_RESP_WOH_COMB))
        resp <- cval_imp$OPTIMUM_RESP_WOH_COMB
      } else if(resp_type == "OPTIMUM_RESP_WH_COMB"){
        cval_imp <- cval_imp %>% filter(!is.na(OPTIMUM_RESP_WH_COMB))
        resp <- cval_imp$OPTIMUM_RESP_WH_COMB
      } else if(resp_type == "OPTIMUM_RESP_WOH_COMB_14d"){
        cval_imp <- cval_imp %>% filter(!is.na(OPTIMUM_RESP_WOH_COMB_14d))
        resp <- cval_imp$OPTIMUM_RESP_WOH_COMB_14d
      } else if(resp_type == "OPTIMUM_RESP_WH_COMB_14d"){
        cval_imp <- cval_imp %>% filter(!is.na(OPTIMUM_RESP_WH_COMB_14d))
        resp <- cval_imp$OPTIMUM_RESP_WH_COMB_14d
      } else{
        stop("No response selected")
      }
      
      a1 <- paste0("glm(resp ~ ")
      a2 <- paste0("", var_to_add_adj, "", collapse = " + ")
      
      ww <- paste0("", a1, "", a2, ", data = cval_imp, 
                     family = binomial(link = 'logit'), maxit = 200)") # remove NAs by default
      fom <- try(eval(parse(text = ww)))
      
      if(method == "GLM + Penalty + Boot"){
        
        # extract model matrix
        mf <- model.matrix(fom)[, -1]
        my <- resp
        
        
        set.seed(AA + 9999) 
        
        cval_imp$missing_wt <- ifelse(cval_imp$missing_wt >  quantile(cval_imp$missing_wt, 0.95, na.rm = T),  quantile(cval_imp$missing_wt, 0.95, na.rm  = T), cval_imp$missing_wt) # clipped at 95% weight value
        
        system.time(cv.fit <- cv.glmnet(mf, my, weights = cval_imp$missing_wt, family = family_type, nfolds = ifelse(is.null(minnum) == FALSE, ceiling(nrow(cval_imp) / minnum), 
                                                                                                                     nfolds), type.measure = measure_type, #"class", #"auc", #"class", 
                                        gamma = seq(0, 1, length.out = ngamma), alpha = alpha_val, relax = relax_me, parallel = FALSE))
        beta.fit[[ll]] <- coef(cv.fit, s = "lambda.min") 
        name_coef <- rownames(beta.fit[[ll]])
        
        
        cores <- detectCores()
        cl <- makeCluster(cores[1] - 10) #not to overload your computer
        registerDoParallel(cl)
        
        system.time(
          betaB <- foreach(bb = 1 : B, .combine = cbind, .packages = c('glmnet',  'ranger', 'doParallel', 'survival')) %dopar% {
            #try(for(bb in 1 : B){
            set.seed(bb + 888 * ll + AA)
            n <- nrow(mf)
            wt <- rexp(n, 1)
            wt <- wt / sum(wt)
          
            system.time(cv.fitB <- cv.glmnet(mf, my, family = family_type, weights = cval_imp$missing_wt * wt, 
                                             nfolds = ifelse(is.null(minnum) == FALSE, ceiling(nrow(cval_imp) / minnum), nfolds), 
                                             type.measure = measure_type, alpha = alpha_val,   
                                             gamma = seq(0, 1, length.out = ngamma), relax = relax_me, parallel = FALSE))
            coef(cv.fitB, s = "lambda.min") 
          })
        
        
        stopCluster(cl)
        rn <- rownames(beta.fit[[ll]])
        
        ## bias-adjusted based on fractional-random weight bootstrap https://par.nsf.gov/servlets/purl/10155761
        boot_bound <- do.call(rbind, lapply(seq_len(length(beta.fit[[ll]])), function(jj) {
          ic <- which(rn == rn[jj]) 
          icb <- which(rownames(betaB) == rn[jj])
          bb <- sort(betaB[icb, ])
          q <- length(which(bb < beta.fit[[ll]][ic])) / B
          lb_ix <- round(B * pnorm((2 * qnorm(q, 0, 1)) + qnorm(alpha / 2, 0, 1)))
          lb_ix <- ifelse(lb_ix == 0, 1, lb_ix)
          ub_ix <- round(B * pnorm((2 * qnorm(q, 0, 1)) + qnorm(1 - (alpha / 2), 0, 1)))
         
          data.frame("names" =  rn[jj], "lb_bias2" = bb[lb_ix], "ub_bias2" = bb[ub_ix])
        }))
        
        
        bias <- (2 * as.numeric(beta.fit[[ll]])) - apply(betaB, 1, mean) 
        bm[[ll]] <- betaB
        
     
        dat_imp[[ll]] <- data.frame("Variable" = name_coef,
                                    "bOR" = exp(as.numeric(beta.fit[[ll]])),
                                    "lb_025" = exp(apply(betaB, 1, finq025)),
                                    "ub_975" = exp(apply(betaB, 1, finq975)),
                                    "sd" = sqrt(apply(betaB, 1, var)),
                                    
                                    "p_val" =  2 * (1 - pnorm(abs(as.numeric(beta.fit[[ll]]) / sqrt(apply(betaB, 1, var))))),
                                    
                                    "lb_bias" = exp(bias - ME * sqrt(apply(betaB, 1, var)) ),
                                    "ub_bias" = exp(bias + ME * sqrt(apply(betaB, 1, var)) ),
                                    
                                    "lb_bias2" = exp(boot_bound$lb_bias2), 
                                    "ub_bias2" = exp(boot_bound$ub_bias2),
                                    
                                    "lb_emp" = exp(2 * as.numeric(beta.fit[[ll]]) - apply(betaB, 1, finq975)),
                                    "ub_emp" = exp(2 * as.numeric(beta.fit[[ll]]) + apply(betaB, 1, finq025)),
                                    
                                    "lb" = exp(as.numeric(beta.fit[[ll]]) - ME * sqrt(apply(betaB, 1, var)) ),  
                                    "ub" = exp(as.numeric(beta.fit[[ll]]) + ME * sqrt(apply(betaB, 1, var)) ))
        
      }
      
    } # ends for each sample
    
    var_to_add <- var_to_add_adj
    
    
    ## bias adjusted based on fractional-random weight bootstrap https://par.nsf.gov/servlets/purl/10155761
    rn <- rownames(beta.fit[[ll]])
    
    gt <- paste0("beta.fit[[", 1 : nsample, "]]", collapse  = " + ")
    gt <- eval(parse(text = gt))
    beta.fit2 <- gt / nsample
    
    dat <- data.frame("Variable" = name_coef,
                      
                      "bOR" = exp(as.numeric(beta.fit2)),
                      
                      "lb_025" = exp(unlist(lapply(seq_len(length(name_coef)), function(jj) { 
                        finq025(unlist(lapply(seq_len(nsample), function(ll) {  bm[[ll]][jj, ] })))
                      }))),
                      
                      "ub_975" = exp(unlist(lapply(seq_len(length(name_coef)), function(jj) { 
                        finq975(unlist(lapply(seq_len(nsample), function(ll) {  bm[[ll]][jj, ] })))
                      }))))
    dat$lb_bias2 <- dat$lb <- dat$lb_025 
    dat$ub_bias2 <- dat$ub <- dat$ub_975 
  }
  
  # survival without imputation
  if(logistic == FALSE & imputeme == FALSE){
    
    if(relax_me == FALSE){
      cval$time <- cval$time # scaling time
    } else{
      cval$time <- (cval$time - min(cval$time)) / (max(cval$time) - min(cval$time)) # min-max standardization
    }
    nasum <- function(xx) {sum(is.na(xx))}
    
    ## Specialized feature sets 
    if(length(which(apply(cval[, colnames(cval) %!in% c(grep("FUP30D_", colnames(cval), value = T), 
                                                        grep("OPTIMUM_", colnames(cval), value = T),
                                                        grep("CD7_PRIM", colnames(cval), value = T),
                                                        grep("CD7_MONO", colnames(cval), value = T),
                                                        grep("CD7_BLASTS", colnames(cval), value = T),
                                                        grep("_COMBO", colnames(cval), value = T),
                                                        grep("ate_wt", colnames(cval), value = T), # weights for subjects
                                                        grep("missing_wt", colnames(cval), value = T),  # missing data weights 
                                                        grep("l_", colnames(cval), value = T))], 2, nasum) >= 1)) >= 1){
      stop("Missing labs. Add unknown label")
    }
    
    
    nm0 <- paste0("", var_to_add, "", collapse = " + " )
    part1 <- paste0("Surv(time, censor) ~ ", nm0, "", collapse = " + ")
    
    ### Survival regression
    dist_select <- ph_check <- NULL
    
    if(method == "Cox-PH + Penalty + Boot"){
      
      fom <- paste0("survival::coxph(", part1, ", ties = 'breslow', iter.max = 1000, 
                  outer.max = 500, robust = FALSE,
                  singular.ok = TRUE, eps = 1e-1,
                  toler.chol = .Machine$double.eps^.15,
                  data = cval)")
      try(fom <- eval(parse(text =  fom)))
      
      
      # extract model matrix
      mf <- mf_back <- model.matrix(fom)
      
      yss <- Surv(cval$time, cval$censor)
      
      cval$missing_wt <- ifelse(cval$missing_wt >  quantile(cval$missing_wt, 0.95, na.rm = T),  quantile(cval$missing_wt, 0.95, na.rm = T), cval$missing_wt) # clipped at 95% weight value to
      
      set.seed(AA) 
      if(psm == TRUE & psm_method == "IPTW"){ # if models are adjusted for propensity scores for treatment effects for a case when there is a treatment variable
        
        cval$ate_wt1 <- ifelse(cval$ate_wt1 >  quantile(cval$ate_wt1, 0.95, na.rm = T),  quantile(cval$ate_wt1, 0.95, na.rm = T), cval$ate_wt1)
        system.time(cv.fit <- cv.glmnet(mf, yss, weights = cval$missing_wt * cval$ate_wt1, 
                                        family = family_type, nfolds = ifelse(is.null(minnum) == FALSE, ceiling(nrow(cval) / minnum), 
                                                                              nfolds), type.measure = measure_type, alpha = alpha_val,   
                                        gamma = seq(0, 1, length.out = ngamma), relax = relax_me, parallel = FALSE))
      } else if(psm == TRUE & psm_method == "splines"){ # if models are adjusted for propensity scores for treatment effects if there is a treatment variable
        mf <- cbind(mf_back, bs_mat)
        
        system.time(cv.fit <- cv.glmnet(mf, yss, weights = cval$missing_wt, 
                                        family = family_type, nfolds = ifelse(is.null(minnum) == FALSE, ceiling(nrow(cval) / minnum), 
                                                                              nfolds), type.measure = measure_type, alpha = alpha_val,   
                                        gamma = seq(0, 1, length.out = ngamma), relax = relax_me, parallel = FALSE))
      } else if(psm == FALSE){
        system.time(cv.fit <- cv.glmnet(mf, yss, weights = cval$missing_wt, 
                                        family = family_type, nfolds = ifelse(is.null(minnum) == FALSE, ceiling(nrow(cval) / minnum), nfolds), type.measure = measure_type, alpha = alpha_val,   
                                        gamma = seq(0, 1, length.out = ngamma), relax = relax_me, parallel = FALSE))
      }  
      
      
      beta.fit <- coef(cv.fit, s = "lambda.min") # This is equivalent ot glmnet(mf, yss, family = "cox", lambda = cv.fit$lambda.min)
      name_coef <- name_coef2 <- rownames(beta.fit)
      name_coef <- name_coef[name_coef2 %!in% grep("spline", name_coef2, value = T)]
      beta.fit <- beta.fit[name_coef2 %!in% grep("spline", name_coef2, value = T), 1]
      time <- cval$time
      censor <- cval$censor 
      t.unique <- sort(unique(time[cval$censor == 1L]))
      centered = FALSE 
      
      if(nnt == TRUE){ ## Marginal risk + NNT calculation by treatment type 
        nntf <- nnt_find(xmat = mf, vecy = yss, var = "TX_TYPE", itm = c(30, 90, 180, 365, 540), wtt = cval$missing_wt)
      }
      
      
      ## Parallel computation for FRWB runs 
      cores <- detectCores()
      cl <- makeCluster(cores[1] - 10) # not to overload computer
      registerDoParallel(cl)
      
      system.time(
        betaB <- foreach(bb = 1 : B, .combine = cbind, .packages = c('glmnet',  'ranger', 'doParallel', 'survival')) %dopar% {
        
          set.seed(bb + 787878 + AA)
          n <- nrow(cval)
          wt <- rexp(n, 1)
          wt <- wt / sum(wt)
          
          if(psm == TRUE & psm_method == "IPTW"){
            system.time(cv.fitB <- cv.glmnet(mf, yss, weights = cval$missing_wt * wt * cval$ate_wt1, 
                                             family = family_type, nfolds = ifelse(is.null(minnum) == FALSE, ceiling(nrow(cval) / minnum), 
                                                                                   nfolds), type.measure = measure_type, alpha = alpha_val,   
                                             gamma = seq(0, 1, length.out = ngamma), relax = relax_me, parallel = FALSE))
            
          } else if(psm == TRUE & psm_method == "splines"){
            mf <- cbind(mf_back, bs_mat)
            system.time(cv.fitB <- cv.glmnet(mf, yss, weights = cval$missing_wt * wt, 
                                             family = family_type, nfolds = ifelse(is.null(minnum) == FALSE, ceiling(nrow(cval) / minnum), 
                                                                                   nfolds), type.measure = measure_type, alpha = alpha_val,   
                                             gamma = seq(0, 1, length.out = ngamma), relax = relax_me, parallel = FALSE))
          } else if(psm == FALSE){  
            system.time(cv.fitB <- cv.glmnet(mf, yss, family = family_type, weights = cval$missing_wt * wt, 
                                             nfolds = ifelse(is.null(minnum) == FALSE, ceiling(nrow(cval) / minnum), nfolds), type.measure = measure_type, alpha = alpha_val,   
                                             gamma = seq(0, 1, length.out = ngamma), relax = relax_me, parallel = FALSE))
          }  
          coef(cv.fitB, s = "lambda.min") 
        })
      
      stopCluster(cl)
      
      
      rn <- name_coef
      betaB <- betaB[name_coef2 %!in% grep("spline", name_coef2, value = T), ]
      ## biase adjusted based on fractional-random weight bootstrap https://par.nsf.gov/servlets/purl/10155761
      boot_bound <- do.call(rbind, lapply(seq_len(length(beta.fit)), function(jj) {
        ic <- which(rn == rn[jj]) 
        icb <- which(rownames(betaB) == rn[jj])
        bb <- sort(betaB[icb, ])
        q <- length(which(bb < beta.fit[ic])) / B
        lb_ix <- round(B * pnorm((2 * qnorm(q, 0, 1)) + qnorm(alpha / 2, 0, 1)))
        lb_ix <- ifelse(lb_ix == 0, 1, lb_ix)
        ub_ix <- round(B * pnorm((2 * qnorm(q, 0, 1)) + qnorm(1 - (alpha / 2), 0, 1)))
        ub_ix <- ifelse(ub_ix == 0, 1, ub_ix)
        data.frame("names" =  rn[jj], "lb_bias2" = bb[lb_ix], "ub_bias2" = bb[ub_ix])
      }))
      
      
      
      bias <- (2 * as.numeric(beta.fit)) - apply(betaB, 1, mean) # following https://math.montana.edu/jobo/thainp/boot.pdf
      
      dat <- data.frame("Variable" = name_coef,
                        "bOR" = exp(as.numeric(beta.fit)),
                        "lb_025" = exp(apply(betaB, 1, finq025)),
                        "ub_975" = exp(apply(betaB, 1, finq975)),
                        "sd" = sqrt(apply(betaB, 1, var)),
                        
                        "p_val" =  2 * (1 - pnorm(abs(as.numeric(beta.fit) / sqrt(apply(betaB, 1, var))))), 
                        
                        "lb_bias" = exp(bias - ME * sqrt(apply(betaB, 1, var)) ),
                        "ub_bias" = exp(bias + ME * sqrt(apply(betaB, 1, var)) ),
                        
                        "lb_bias2" = exp(boot_bound$lb_bias2), 
                        "ub_bias2" = exp(boot_bound$ub_bias2),
                        
                        "lb_emp" = exp(2 * as.numeric(beta.fit) - apply(betaB, 1, finq975)),
                        "ub_emp" = exp(2 * as.numeric(beta.fit) + apply(betaB, 1, finq025)),
                        "lb" = exp(as.numeric(beta.fit) - ME * sqrt(apply(betaB, 1, var)) ),  
                        "ub" = exp(as.numeric(beta.fit) + ME * sqrt(apply(betaB, 1, var)) ))
      
    } 
 
  }
  
  # survival with imputation
  if(logistic == FALSE & imputeme == TRUE){
    
    dat_imp <- bm <- beta.fit <- trt.fitb <- rmst.fitb <- varimpb <- trt.fit <- rmst.fit <- varimp  <- list(); 
    
    dat_orig <- dat <- trt.fitb <- list()
    
    for(ll in 1 : nsample){ # run for each imputed dataset 
      
      
      var_to_add_adj <- var_to_add[var_to_add %!in% mc$dropme]
      cval_imp <- mc$subd_imp[[ll]][, colnames(mc$subd_imp[[ll]]) %!in% mc$dropme]
      
      if(relax_me == FALSE){
        cval_imp$time <- cval_imp$time # scaling time
      } else{
        cval_imp$time <- (cval_imp$time - min(cval_imp$time)) / (max(cval_imp$time) - min(cval_imp$time)) # min-max standardization
      }
      
      nasum <- function(xx) {sum(is.na(xx))}
      
      if(length(which(apply(cval_imp[, colnames(cval_imp) %!in% c(grep("FUP30D_", colnames(cval_imp), value = T), 
                                                                  grep("OPTIMUM_", colnames(cval_imp), value = T),
                                                                  grep("CD7_PRIM", colnames(cval), value = T),
                                                                  grep("CD7_MONO", colnames(cval), value = T),
                                                                  grep("CD7_BLASTS", colnames(cval), value = T),
                                                                  grep("_COMBO", colnames(cval), value = T),
                                                                  grep("ate_wt", colnames(cval_imp), value = T),
                                                                  grep("l_", colnames(cval_imp), value = T))], 2, nasum) >= 1)) >= 1){
        stop("Missing labs. Add unknown label")
      }
      
      if(psm == TRUE & psm_method == "splines"){
        library(splines2)
        
        try(l_vc <- cval_imp$l_vc)
        if(sum(is.na(l_vc)) != length(l_vc)){
          if(length(which(sort(l_vc) == max(sort(l_vc)))) > 1 | length(which(sort(l_vc) == min(sort(l_vc)))) > 1){ # adjustment to by-pass numerical error at the boundary
            try(knots1 <- seq(head(sort(l_vc), 2)[2], tail(sort(l_vc), 2)[1] - 0.001, length.out = 6))  
          } else{
            try(knots1 <- seq(head(sort(l_vc), 2)[2], tail(sort(l_vc), 2)[1], length.out = 6))
          }
          if(min(knots1) <= min(l_vc)){
            knots1[1] <- knots1[1] + 0.001
          }
          if(max(knots1) <= max(l_vc)){
            knots1[length(knots1)] <- knots1[length(knots1)] - 0.001
          }
          try(bs_l_vc <- bSpline(l_vc, knots = knots1, degree = 0, intercept = TRUE))
        }
        
        try(l_vac <- cval_imp$l_vac)
        if(sum(is.na(l_vac)) != length(l_vac)){
          if(length(which(sort(l_vac) == max(sort(l_vac)))) > 1 | length(which(sort(l_vac) == min(sort(l_vac)))) > 1){ # adustment to by-pass error
            try(knots2 <- seq(head(sort(l_vac), 2)[2], tail(sort(l_vac), 2)[1] - 0.001, length.out = 6))  
          } else{
            try(knots2 <- seq(head(sort(l_vac), 2)[2], tail(sort(l_vac), 2)[1], length.out = 6))
          }
          if(min(knots2) <= min(l_vac)){
            knots2[1] <- knots2[1] + 0.001
          }
          if(max(knots2) <= max(l_vac)){
            knots2[length(knots2)] <- knots2[length(knots2)] - 0.001
          }
          
          try(bs_l_vac <- bSpline(l_vac, knots = knots2, degree = 0, intercept = TRUE))
        }
        
        try(l_a75c10 <- cval_imp$l_a75c10)
        if(sum(is.na(l_a75c10)) != length(l_a75c10)){
          if(length(which(sort(l_a75c10) == max(sort(l_a75c10)))) > 1 | length(which(sort(l_a75c10) == min(sort(l_a75c10)))) > 1){ # adustment to by-pass error
            try(knots3 <- seq(head(sort(l_a75c10), 2)[2], tail(sort(l_a75c10), 2)[1] - 0.001, length.out = 6))  
          } else{
            try(knots3 <- seq(head(sort(l_a75c10), 2)[2], tail(sort(l_a75c10), 2)[1], length.out = 6))
          }
          if(min(knots3) <= min(l_a75c10)){
            knots3[1] <- knots3[1] + 0.001
          }
          if(max(knots3) <= max(l_a75c10)){
            knots3[length(knots3)] <- knots3[length(knots3)] - 0.001
          }
          try(bs_l_a75c10 <- bSpline(l_a75c10, knots = knots3, degree = 0, intercept = TRUE))
        }
        
        try(l_a75c20 <- cval_imp$l_a75c20)
        if(sum(is.na(l_a75c20)) != length(l_a75c20)){
          if(length(which(sort(l_a75c20) == max(sort(l_a75c20)))) > 1 | length(which(sort(l_a75c20) == min(sort(l_a75c20)))) > 1){ # adustment to by-pass error
            try(knots4 <- seq(head(sort(l_a75c20), 2)[2], tail(sort(l_a75c20), 2)[1] - 0.001, length.out = 6))  
          } else{
            try(knots4 <- seq(head(sort(l_a75c20), 2)[2], tail(sort(l_a75c20), 2)[1], length.out = 6))
          }
          if(min(knots4) <= min(l_a75c20)){
            knots4[1] <- knots4[1] + 0.001
          }
          if(max(knots4) <= max(l_a75c20)){
            knots4[length(knots4)] <- knots4[length(knots4)] - 0.001
          }
          try(bs_l_a75c20 <- bSpline(l_a75c20, knots = knots4, degree = 0, intercept = TRUE))
        }
        
        
        if(alltx == "all") {#%in% c("ven/aza", "ven/aza/cusa", "aza75/cusa20", "aza75/cusa10", "ven/cusa") == TRUE) )
          bs_mat <- cbind(bs_l_vc, bs_l_vac, bs_l_a75c10, bs_l_a75c20)
        } else if(alltx == "a75c20 vs va"){ # tx-1 vs tx-2
          bs_mat <- cbind(bs_l_a75c20)
        } else if(alltx == "a75c10 vs va"){ # tx-3 vs tx-2
          bs_mat <- cbind(bs_l_a75c10)
        }
        
        colnames(bs_mat) <- paste0("spline", 1 : ncol(bs_mat), "")
      }
      
      nm0 <- paste0("", var_to_add_adj, "", collapse = " + " )
      part1 <- paste0("Surv(time, censor) ~ ", nm0, "", collapse = " + ")
      
      if(method == "Cox-PH + Penalty + Boot"){
        
        fom <- paste0("survival::coxph(", part1, ", ties = 'breslow', iter.max = 350, 
                  outer.max = 500, robust = TRUE,
                  singular.ok = TRUE, eps = 1e-1,
                  toler.chol = .Machine$double.eps^.35,
                  data = cval_imp)")
        try(fom <- eval(parse(text =  fom)))
        
        # extract model matrix
        mf <- mf_back <- model.matrix(fom)
        
        yss <- Surv(cval_imp$time, cval_imp$censor)
        
        cval_imp$missing_wt <- ifelse(cval_imp$missing_wt >  quantile(cval_imp$missing_wt, 0.95, na.rm = T),  quantile(cval_imp$missing_wt, 0.95, na.rm  = T), cval_imp$missing_wt) # clipped at 95% weight value
        
        set.seed(AA + 111 * ll) 
        if(psm == TRUE & psm_method == "IPTW"){
          cval_imp$ate_wt1 <- ifelse(cval_imp$ate_wt1 >  quantile(cval_imp$ate_wt1, 0.95, na.rm = T),  quantile(cval_imp$ate_wt1, 0.95, na.rm = T), cval_imp$ate_wt1)
          
          system.time(cv.fit <- cv.glmnet(mf, yss, weights = cval_imp$missing_wt * cval_imp$ate_wt1, 
                                          family = family_type, nfolds = ifelse(is.null(minnum) == FALSE, ceiling(nrow(cval_imp) / minnum), 
                                                                                nfolds), type.measure = measure_type, alpha = alpha_val,   
                                          gamma = seq(0, 1, length.out = ngamma), relax = relax_me, parallel = FALSE))
        } else if(psm == TRUE & psm_method == "OW"){
          
          system.time(cv.fit <- cv.glmnet(mf, yss, weights = cval_imp$missing_wt * cval_imp$ate_ow, 
                                          family = family_type, nfolds = ifelse(is.null(minnum) == FALSE, ceiling(nrow(cval_imp) / minnum), 
                                                                                nfolds), type.measure = measure_type, alpha = alpha_val,   
                                          gamma = seq(0, 1, length.out = ngamma), relax = relax_me, parallel = FALSE))
          
        } else if(psm == TRUE & psm_method == "SMR"){
          
          cval_imp$ate_smr <- ifelse(cval_imp$ate_smr >  quantile(cval_imp$ate_smr, 0.85, na.rm = T),  quantile(cval_imp$ate_smr, 0.85, na.rm = T), cval_imp$ate_smr)
          
          system.time(cv.fit <- cv.glmnet(mf, yss, weights = cval_imp$missing_wt * cval_imp$ate_smr, 
                                          family = family_type, nfolds = ifelse(is.null(minnum) == FALSE, ceiling(nrow(cval_imp) / minnum), 
                                                                                nfolds), type.measure = measure_type, alpha = alpha_val,   
                                          gamma = seq(0, 1, length.out = ngamma), relax = relax_me, parallel = FALSE))
          
        } else if(psm == TRUE & psm_method == "splines"){
          mf <- cbind(mf_back, bs_mat)
          system.time(cv.fit <- cv.glmnet(mf, yss, weights = cval_imp$missing_wt, 
                                          family = family_type, nfolds = ifelse(is.null(minnum) == FALSE, ceiling(nrow(cval_imp) / minnum), 
                                                                                nfolds), type.measure = measure_type, alpha = alpha_val,   
                                          gamma = seq(0, 1, length.out = ngamma), relax = relax_me, parallel = FALSE))
        } else if(psm == FALSE){
          system.time(cv.fit <- cv.glmnet(mf, yss, weights = cval_imp$missing_wt, 
                                          family = family_type, nfolds = ifelse(is.null(minnum) == FALSE, ceiling(nrow(cval_imp) / minnum), 
                                                                                nfolds), type.measure = measure_type, alpha = alpha_val,   
                                          gamma = seq(0, 1, length.out = ngamma), relax = relax_me, parallel = FALSE))
        }  
        
        beta.fit[[ll]] <- coef(cv.fit, s = "lambda.min") # This is equivalent ot glmnet(mf, yss, family = "cox", lambda = cv.fit$lambda.min)
        name_coef <- name_coef2 <- rownames(beta.fit[[ll]])
        name_coef <- name_coef[name_coef2 %!in% grep("spline", name_coef2, value = T)]
        beta.fit[[ll]] <- beta.fit[[ll]][name_coef2 %!in% grep("spline", name_coef2, value = T), 1]
        time <- cval_imp$time
        censor <- cval_imp$censor 
        t.unique <- sort(unique(time[cval$censor == 1L]))
        centered = FALSE 
        
        cores <- detectCores()
        cl <- makeCluster(cores[1] - 2) 
        registerDoParallel(cl)
        
        
        system.time(
          betaB <- foreach(bb = 1 : B, .combine = cbind, .packages = c('glmnet',  'ranger', 'doParallel', 'survival')) %dopar% {
           
            set.seed(bb + AA + 111 * ll)
            n <- nrow(cval_imp)
            wt <- rexp(n, 1)
            wt <- wt / sum(wt)
           
            
            if(psm == TRUE & psm_method == "IPTW"){
              system.time(cv.fitB <- cv.glmnet(mf, yss, weights = cval_imp$missing_wt * wt * cval_imp$ate_wt1, 
                                               family = family_type, nfolds = ifelse(is.null(minnum) == FALSE, ceiling(nrow(cval_imp) / minnum), 
                                                                                     nfolds), type.measure = measure_type, alpha = alpha_val,   
                                               gamma = seq(0, 1, length.out = ngamma), relax = relax_me, parallel = FALSE))
              
            } else if(psm == TRUE & psm_method == "OW"){
              system.time(cv.fitB <- cv.glmnet(mf, yss, weights = cval_imp$missing_wt * wt * cval_imp$ate_ow, 
                                               family = family_type, nfolds = ifelse(is.null(minnum) == FALSE, ceiling(nrow(cval_imp) / minnum), 
                                                                                     nfolds), type.measure = measure_type, alpha = alpha_val,   
                                               gamma = seq(0, 1, length.out = ngamma), relax = relax_me, parallel = FALSE))
              
            } else if(psm == TRUE & psm_method == "SMR"){
              system.time(cv.fitB <- cv.glmnet(mf, yss, weights = cval_imp$missing_wt * wt * cval_imp$ate_smr, 
                                               family = family_type, nfolds = ifelse(is.null(minnum) == FALSE, ceiling(nrow(cval_imp) / minnum), 
                                                                                     nfolds), type.measure = measure_type, alpha = alpha_val,   
                                               gamma = seq(0, 1, length.out = ngamma), relax = relax_me, parallel = FALSE))
              
            } else if(psm == TRUE & psm_method == "splines"){
              mf <- cbind(mf_back, bs_mat)
              system.time(cv.fitB <- cv.glmnet(mf, yss, weights = cval_imp$missing_wt * wt, 
                                               family = family_type, nfolds = ifelse(is.null(minnum) == FALSE, ceiling(nrow(cval_imp) / minnum), 
                                                                                     nfolds), type.measure = measure_type, alpha = alpha_val,   
                                               gamma = seq(0, 1, length.out = ngamma), relax = relax_me, parallel = FALSE))
            } else if(psm == FALSE){  
              system.time(cv.fitB <- cv.glmnet(mf, yss, family = family_type, weights = cval_imp$missing_wt * wt, 
                                               nfolds = ifelse(is.null(minnum) == FALSE, ceiling(nrow(cval_imp) / minnum), nfolds), type.measure = measure_type, alpha = alpha_val,   
                                               gamma = seq(0, 1, length.out = ngamma), relax = relax_me, parallel = FALSE))
            }  
            
            coef(cv.fitB, s = "lambda.min") 
          })
        
        stopCluster(cl)
        
        
        rn <- name_coef
        betaB <- betaB[name_coef2 %!in% grep("spline", name_coef2, value = T), ]
        ## bias-adjusted based on fractional-random weight bootstrap https://par.nsf.gov/servlets/purl/10155761
        boot_bound <- do.call(rbind, lapply(seq_len(length(beta.fit[[ll]])), function(jj) {
          ic <- which(rn == rn[jj]) 
          icb <- which(rownames(betaB) == rn[jj])
          bb <- sort(betaB[icb, ])
          q <- length(which(bb < beta.fit[[ll]][ic])) / B
          lb_ix <- round(B * pnorm((2 * qnorm(q, 0, 1)) + qnorm(alpha / 2, 0, 1)))
          lb_ix <- ifelse(lb_ix == 0, 1, lb_ix)
          ub_ix <- round(B * pnorm((2 * qnorm(q, 0, 1)) + qnorm(1 - (alpha / 2), 0, 1)))
          #c( bb[lb_ix],  bb[ub_ix])
          data.frame("names" =  rn[jj], "lb_bias2" = bb[lb_ix], "ub_bias2" = bb[ub_ix])
        }))
        
        
        
        bias <- (2 * as.numeric(beta.fit[[ll]])) - apply(betaB, 1, mean) # following https://math.montana.edu/jobo/thainp/boot.pdf
       
        bm[[ll]] <- betaB
        
        
        dat_imp[[ll]] <- data.frame("Variable" = name_coef,
                                    "bOR" = exp(as.numeric(beta.fit[[ll]])),
                                    "lb_025" = exp(apply(betaB, 1, finq025)),
                                    "ub_975" = exp(apply(betaB, 1, finq975)),
                                    "sd" = sqrt(apply(betaB, 1, var)),
                                    "p_val" =  2 * (1 - pnorm(abs(as.numeric(beta.fit[[ll]]) / sqrt(apply(betaB, 1, var))))),
                                    "lb_bias" = exp(bias - ME * sqrt(apply(betaB, 1, var)) ),
                                    "ub_bias" = exp(bias + ME * sqrt(apply(betaB, 1, var)) ),
                                    
                                    "lb_bias2" = exp(boot_bound$lb_bias2), 
                                    "ub_bias2" = exp(boot_bound$ub_bias2),
                                    
                                    "lb_emp" = exp(2 * as.numeric(beta.fit[[ll]]) - apply(betaB, 1, finq975)),
                                    "ub_emp" = exp(2 * as.numeric(beta.fit[[ll]]) + apply(betaB, 1, finq025)),
                                   
                                    "lb" = exp(as.numeric(beta.fit[[ll]]) - ME * sqrt(apply(betaB, 1, var)) ),  
                                    "ub" = exp(as.numeric(beta.fit[[ll]]) + ME * sqrt(apply(betaB, 1, var)) ))
        
      } # methods ends 
      
    }  # for each imp dataset 
    
    var_to_add <- var_to_add_adj
    
    ## biase based on fractional-random weight bootstrap https://par.nsf.gov/servlets/purl/10155761
    if(method %!in% c("CausalML", "CausalML1000", "CausalML120", "CausalML365", grep("Causal", method, value = T))){
      rn <- rownames(beta.fit[[ll]])
      
      gt <- paste0("beta.fit[[", 1 : nsample, "]]", collapse  = " + ")
      gt <- eval(parse(text = gt))
      beta.fit2 <- gt / nsample
      
      dat <- data.frame("Variable" = name_coef,
                        
                        "bOR" = exp(as.numeric(beta.fit2)),
                        
                        "lb_025" = exp(unlist(lapply(seq_len(length(name_coef)), function(jj) { 
                          finq025(unlist(lapply(seq_len(nsample), function(ll) {  bm[[ll]][jj, ] })))
                        }))),
                        
                        "ub_975" = exp(unlist(lapply(seq_len(length(name_coef)), function(jj) { 
                          finq975(unlist(lapply(seq_len(nsample), function(ll) {  bm[[ll]][jj, ] })))
                        }))))
      dat$lb_bias2 <- dat$lb <- dat$lb_025 
      dat$ub_bias2 <- dat$ub <- dat$ub_975 
    }
    
  } # loop closes
  
  
  dl <- NULL
  
  if(imputeme ==  FALSE){
    cval <- cval %>% mutate_if(is.factor, as.character)
  }
  
  if(imputeme == TRUE){
    cval <- cval_imp
    cval <- cval %>% mutate_if(is.character, as.factor)
  }
  
  
  if(delete_NP_param == TRUE){
    ## Start here deleting NP
    dl <- unlist(lapply(seq_len(nrow(dat)), function(ii) {
      length(grep("Missing", as.character(dat$Variable[ii]), value = F))}))
    ixx <- which(dl >= 1)
    if(length(ixx) > 0){
      dat <- dat[-ixx, ]
    }
    dat <- dat %>% mutate_if(is.factor, as.character)
    dat <- dat %>% mutate_if(is.character, as.factor)
  }
  
  
  output <- list (dat = dat, 
                  var_to_add = var_to_add, 
                  data = cval)

}



## fitting counterfactual estimations
# Inputs are as same as depicted in Inter_TX_multivariate_profile_paperML() except 
# nnt_var = Variables for which marginal risks for (+ve) and (-ve) need to be computed
# nnt_tm = Time at which marginal risks were evaluated (a subset of grid values in dynamic profiles) 
Inter_TX_multivariate_profile_counterfact_paperML <- function(udata, tx_type, method, var_to_add,
                                                              what = "N", 
                                                              alpha = 0.05,
                                                              delete_np = FALSE,
                                                              delete_NP_param = TRUE,
                                                              logistic = FALSE,
                                                              
                                                              nnt = TRUE,
                                                              nnt_var = c("NGS_KRAS", "NGS_NRAS"), 
                                                              nnt_tm = c(30, 90, 180, 365, 540),
                                                              
                                                              psm = FALSE, 
                                                              psm_method = "IPTW", 
                                                              
                                                              psm_model= NULL, 
                                                              psm_var = c("DEM_AGE", "ECOG"),
                                                              alltx = "drop vc", 
                                                              
                                                              response = 30,
                                                              resp_type = "OS",
                                                              B = 100,
                                                              nfolds = 15, 
                                                              ngamma = 30, 
                                                              relax_me = FALSE, 
                                                              AA = 1200,
                                                              minnum = 5,
                                                              alpha_val = 0,
                                                              measure_type = "deviance",
                                                              family_type = "cox",
                                                              imputeme = FALSE,
                                                              nsample = 5) {
  
  udata <- as.data.frame(udata)
  
  
  if(is.list(var_to_add) == TRUE){
    var_to_impute <- var_to_add[[2]]
    var_to_add <- var_to_add[[1]]
  } else{
    var_to_impute <- var_to_add
  }
  
  
  symbol <- "\u2265"; symbol2 <- intToUtf8(8804)
  ME = abs(qnorm(alpha / 2))
  
  ## deleting NPs and its associated subject for better interpretation
  udata <- udata %>% mutate_if(is.factor, as.character) # deleting original factoring labels
  
  if(delete_np == TRUE){
    udata <- udata %>% na_if("NP") %>% na_if("Not performed") %>% na_if("Unable to assess") %>% na_if("Unknown")  %>% na_if("Missing")
    ## deleting patients for NPs
    adj2 <- which(as.vector(apply(udata, 1, function(ff) sum(is.na(ff)))) >= 1) 
    if(length(adj2) >= 1){
      udata <- udata[- adj2, ]
    }
  }
  
  ## Categorizing numerical lab variables
  if(length(grep("RANGE", var_to_add, value = T)) > 0){ 
    try(cval <- numeric_reference_profile_ov(udata)$result)
  }
  
  
  finq025 <- function(ff) {quantile(ff, 0.025)}
  finq975 <- function(ff) {quantile(ff, 0.975)}
  finq050 <- function(ff) {quantile(ff, 0.5)}
  
  
  
  
  cval <- cval[, colnames(cval) %in% c("MRN", "time", "censor", "TX_TYPE", 
                                       "FUP30D_RESP", 
                                       "OPTIMUM_RESP_WOH_COMB", "OPTIMUM_RESP_WH_COMB",
                                       "OPTIMUM_RESP_WOH_COMB_14d", "OPTIMUM_RESP_WH_COMB_14d", 
                                       var_to_add, var_to_impute, resp_type, psm_var, "missing_wt", "ate_wt1")]
  
  
  
  ## convert factor into character
  cval <- cval %>% mutate_if(is.character, as.factor) # re-creating factors
  try(cval <-  categorical_reference_profile_ov(udata = cval)$result)
  
  
  ## If imputed models are needed & whether propensity score models are needed to be fitted 
  if(imputeme == TRUE & psm == FALSE){  
    system.time(mc <- multivariate_impute_ov(udata, varb = var_to_impute, delete_np, resp_type, nsample = nsample, AA = AA, psm = FALSE, psm_var = psm_var, alltx = alltx, psm_model= psm_model))
  } else if(imputeme == TRUE & psm == TRUE & is.null(psm_model) == FALSE){ # Generated psm weights within MICE step for each imputed datatset 
    system.time(mc <- multivariate_impute_ov(udata, varb = var_to_impute, delete_np, resp_type, nsample = nsample, AA = AA, psm = TRUE, psm_var = psm_var, alltx = alltx, psm_model= psm_model))
  }
  
  
  ## checking if each category (except missing) label in var_to_add has at-least minn values 
  notry <- grep("_COMBO", var_to_add, value = T)
  dll <- which(
    unlist(lapply(seq_len(length(var_to_add)), function(jj) {
      
      QQ <- paste0("table(cval$", var_to_add[jj], ")")
      QQ <- eval(parse(text = QQ))
      
      if(var_to_add[jj] %!in% notry & length(which(QQ[names(QQ) %!in% "Missing"] < minnum)) >= 1){
        1111 
      } else{
        0
      }
      
    })) == 1111)
  if(length(dll) >= 1){
    nnt_var <- nnt_var[ nnt_var %!in% var_to_add[dll] ]
    var_to_add <- var_to_add[- dll] # delete these and update var_to_add 
    var_to_impute <- var_to_impute[- dll]
    
  }
  
  # survival without imputation
  if(logistic == FALSE & imputeme == FALSE){
    
    if(relax_me == FALSE){
      cval$time <- cval$time # scaling time
    } else{
      cval$time <- (cval$time - min(cval$time)) / (max(cval$time) - min(cval$time)) # min-max standardization
    }
    
    nasum <- function(xx) {sum(is.na(xx))}
    
    if(length(which(apply(cval[, colnames(cval) %!in% c(grep("FUP30D_", colnames(cval), value = T), 
                                                        grep("OPTIMUM_", colnames(cval), value = T),
                                                        grep("_COMBO", colnames(cval), value = T),
                                                        grep("ate_wt", colnames(cval), value = T),
                                                        grep("missing_wt", colnames(cval), value = T), 
                                                        grep("l_", colnames(cval), value = T))], 2, nasum) >= 1)) >= 1){
      stop("Missing labs. Add unknown label")
    }
    
    if(psm == TRUE & psm_method == "splines"){
      library(splines2)
      
      try(l_vc <- cval$l_vc)
      if(sum(is.na(l_vc)) != length(l_vc)){
        if(length(which(sort(l_vc) == max(sort(l_vc)))) > 1 | length(which(sort(l_vc) == min(sort(l_vc)))) > 1){ # adjustment to by-pass numerical error
          try(knots1 <- seq(head(sort(l_vc), 2)[2], tail(sort(l_vc), 2)[1] - 0.001, length.out = 6))  
        } else{
          try(knots1 <- seq(head(sort(l_vc), 2)[2], tail(sort(l_vc), 2)[1], length.out = 6))
        }
        if(min(knots1) <= min(l_vc)){
          knots1[1] <- knots1[1] + 0.001
        }
        if(max(knots1) <= max(l_vc)){
          knots1[length(knots1)] <- knots1[length(knots1)] - 0.001
        }
        try(bs_l_vc <- bSpline(l_vc, knots = knots1, degree = 0, intercept = TRUE))
      }
      
      
      try(l_vac <- cval$l_vac)
      if(sum(is.na(l_vac)) != length(l_vac)){
        if(length(which(sort(l_vac) == max(sort(l_vac)))) > 1 | length(which(sort(l_vac) == min(sort(l_vac)))) > 1){ # adjustment to by-pass numerical error
          try(knots2 <- seq(head(sort(l_vac), 2)[2], tail(sort(l_vac), 2)[1] - 0.001, length.out = 6))  
        } else{
          try(knots2 <- seq(head(sort(l_vac), 2)[2], tail(sort(l_vac), 2)[1], length.out = 6))
        }
        if(min(knots2) <= min(l_vac)){
          knots2[1] <- knots2[1] + 0.001
        }
        if(max(knots2) <= max(l_vac)){
          knots2[length(knots2)] <- knots2[length(knots2)] - 0.001
        }
        try(bs_l_vac <- bSpline(l_vac, knots = knots2, degree = 0, intercept = TRUE))
      }
      
      try(l_a75c10 <- cval$l_a75c10)
      if(sum(is.na(l_a75c10)) != length(l_a75c10)){
        if(length(which(sort(l_a75c10) == max(sort(l_a75c10)))) > 1 | length(which(sort(l_a75c10) == min(sort(l_a75c10)))) > 1){ # adjustment to by-pass numerical error
          try(knots3 <- seq(head(sort(l_a75c10), 2)[2], tail(sort(l_a75c10), 2)[1] - 0.001, length.out = 6))  
        } else{
          try(knots3 <- seq(head(sort(l_a75c10), 2)[2], tail(sort(l_a75c10), 2)[1], length.out = 6))
        }
        if(min(knots3) <= min(l_a75c10)){
          knots3[1] <- knots3[1] + 0.001
        }
        if(max(knots3) <= max(l_a75c10)){
          knots3[length(knots3)] <- knots3[length(knots3)] - 0.001
        }
        try(bs_l_a75c10 <- bSpline(l_a75c10, knots = knots3, degree = 0, intercept = TRUE))
      }
      
      try(l_a75c20 <- cval$l_a75c20)
      if(sum(is.na(l_a75c20)) != length(l_a75c20)){
        if(length(which(sort(l_a75c20) == max(sort(l_a75c20)))) > 1 | length(which(sort(l_a75c20) == min(sort(l_a75c20)))) > 1){ # adjustment to by-pass numerical error
          try(knots4 <- seq(head(sort(l_a75c20), 2)[2] + 0.00001, tail(sort(l_a75c20), 2)[1] - 0.001, length.out = 6))  
        } else{
          try(knots4 <- seq(head(sort(l_a75c20), 2)[2], tail(sort(l_a75c20), 2)[1], length.out = 6))
        }
        if(min(knots4) <= min(l_a75c20)){
          knots4[1] <- knots4[1] + 0.001
        }
        if(max(knots4) <= max(l_a75c20)){
          knots4[length(knots4)] <- knots4[length(knots4)] - 0.001
        }
        
        try(bs_l_a75c20 <- bSpline(l_a75c20, knots = knots4, degree = 0, intercept = TRUE))
      }
      
      
      if(alltx == "all") {
        bs_mat <- cbind(bs_l_vc, bs_l_vac, bs_l_a75c10, bs_l_a75c20)
        
      } else if(alltx == "a75c20 vs va"){ # tx-2 vs tx-1
        bs_mat <- cbind(bs_l_a75c20)
      } else if(alltx == "a75c10 vs va"){ # tx-3 vs tx-1
        bs_mat <- cbind(bs_l_a75c10)
      }
      
      
      colnames(bs_mat) <- paste0("spline", 1 : ncol(bs_mat), "")
    }
    
    nm0 <- paste0("", var_to_add, "", collapse = " + " )
    part1 <- paste0("Surv(time, censor) ~ ", nm0, "", collapse = " + ")
    
    
    if(method == "Cox-PH + Penalty + Boot"){
      
      fom <- paste0("survival::coxph(", part1, ", ties = 'breslow', iter.max = 1000, 
                  outer.max = 500, robust = FALSE,
                  singular.ok = TRUE, eps = 1e-1,
                  toler.chol = .Machine$double.eps^.25,
                  data = cval)")
      try(fom <- eval(parse(text =  fom)))
      
      # extract model matrix
      mf <- mf_back <- model.matrix(fom)
      
      
      yss <- Surv(cval$time, cval$censor)
      
      set.seed(AA) 
      if(psm == TRUE & psm_method == "IPTW"){ # if models are adjusted for propensity scores for treatment effects for a case when there is a treatment variable
        
        cval$ate_wt1 <- ifelse(cval$ate_wt1 >  quantile(cval$ate_wt1, 0.95, na.rm = T),  quantile(cval$ate_wt1, 0.95, na.rm = T), cval$ate_wt1)
        system.time(cv.fit <- cv.glmnet(mf, yss, weights = cval$missing_wt * cval$ate_wt1, 
                                        family = family_type, nfolds = ifelse(is.null(minnum) == FALSE, ceiling(nrow(cval) / minnum), 
                                                                              nfolds), type.measure = measure_type, alpha = alpha_val,   
                                        gamma = seq(0, 1, length.out = ngamma), relax = relax_me, parallel = FALSE))
      } else if(psm == TRUE & psm_method == "splines"){ # if models are adjusted for propensity scores for treatment effects if there is a treatment variable
        mf <- cbind(mf_back, bs_mat)
        
        system.time(cv.fit <- cv.glmnet(mf, yss, weights = cval$missing_wt, 
                                        family = family_type, nfolds = ifelse(is.null(minnum) == FALSE, ceiling(nrow(cval) / minnum), 
                                                                              nfolds), type.measure = measure_type, alpha = alpha_val,   
                                        gamma = seq(0, 1, length.out = ngamma), relax = relax_me, parallel = FALSE))
      } else if(psm == FALSE){
        system.time(cv.fit <- cv.glmnet(mf, yss, weights = cval$missing_wt, 
                                        family = family_type, nfolds = ifelse(is.null(minnum) == FALSE, ceiling(nrow(cval) / minnum), nfolds), type.measure = measure_type, alpha = alpha_val,   
                                        gamma = seq(0, 1, length.out = ngamma), relax = relax_me, parallel = FALSE))
      }  
      
      
      beta.fit <- coef(cv.fit, s = "lambda.min") # This is equivalent ot glmnet(mf, yss, family = "cox", lambda = cv.fit$lambda.min)
      name_coef <- name_coef2 <- rownames(beta.fit)
      name_coef <- name_coef[name_coef2 %!in% grep("spline", name_coef2, value = T)]
      beta.fit <- beta.fit[name_coef2 %!in% grep("spline", name_coef2, value = T), 1]
      time <- cval$time
      censor <- cval$censor 
      t.unique <- sort(unique(time[cval$censor == 1L]))
      centered = FALSE 
      
      tm <- inx <- m1 <- m2 <- m3 <- m4 <- m5 <- p1 <- p2 <- p3 <- p4 <- p5 <- NULL
      if(nnt == TRUE){ 
        nntf <- nnt_find(fitb = cv.fit, xmat = mf, vecy = yss, var = nnt_var, itm = nnt_tm, wtt = cval$missing_wt)
        m1 <- nntf$m1; m2 <- nntf$m2; m3 <- nntf$m3; m4 <- nntf$m4; m5 <- nntf$m5; 
        p1 <- nntf$p1; p2 <- nntf$p2; p3 <- nntf$p3; p4 <- nntf$p4; p5 <- nntf$p5; 
        inx <- nntf$inx; 
        try(tm <- nntf$tm)
      }
      
      
      # parallel computation for bootstraps
      comb <- function(x, ...) {
        lapply(seq_along(x),
               function(i) c(x[[i]], lapply(list(...), function(y) y[[i]])))
      }
      
      
      cores <- detectCores()
      cl <- makeCluster(cores[1] - 10) 
      registerDoParallel(cl)
      
      
      system.time(
        rall <- foreach(bb = 1 : B, .combine = 'comb', .multicombine=TRUE, .errorhandling = 'pass',
                        
                        .init = list(list(), list(), list(), list(), list(), list()), 
                        .export = c("nnt_find",  '%!in%'), 
                        .packages = c('glmnet',  'ranger', 'doParallel', 'survival')) %dopar% {
                          
                          unregister <- function() {
                            env <- foreach:::.foreachGlobals
                            rm(list=ls(name=env), pos=env)
                          }    
                          
                          unregister()
                          print(paste0("now-", bb, ""))     
                          
                          
                          #try(for(bb in 1 : B){
                          set.seed(bb + 7878 + AA)
                          n <- nrow(cval)
                          wt <- rexp(n, 1)
                          wt <- wt / sum(wt)
                          #i1 <- sample(seq_len(nn), nn, replace = TRUE)
                          #mfB <- mf[i1, ];  yssB <- Surv(cval$time[i1], cval$censor[i1])
                          #mfB <- mf;  yssB <- Surv(cval$time, cval$censor)
                          if(psm == TRUE & psm_method == "IPTW" & nnt == TRUE){
                            system.time(cv.fitB <- cv.glmnet(mf, yss, weights = cval$missing_wt * wt * cval$ate_wt1, 
                                                             family = family_type, nfolds = ifelse(is.null(minnum) == FALSE, ceiling(nrow(cval) / minnum), 
                                                                                                   nfolds), type.measure = measure_type, alpha = alpha_val,   
                                                             gamma = seq(0, 1, length.out = ngamma), relax = relax_me, parallel = FALSE))
                            
                            
                            nntf <- nnt_find(fitb = cv.fitB, xmat = mf, vecy = yss, var = nnt_var, itm = nnt_tm, wtt = cval$missing_wt * wt * cval$ate_wt1)
                            m1b <- nntf$m1; m2b <- nntf$m2; m3b <- nntf$m3; m4b <- nntf$m4; m5b <- nntf$m5; 
                            p1b <- nntf$p1; p2b <- nntf$p2; p3b <- nntf$p3; p4b <- nntf$p4; p5b <- nntf$p5; 
                            inxb <- nntf$inx; tmb <- nntf$tm
                            
                            
                          } else if(psm == TRUE & psm_method == "splines" & nnt ==  TRUE){
                            mf <- cbind(mf_back, bs_mat)
                            system.time(cv.fitB <- cv.glmnet(mf, yss, weights = cval$missing_wt * wt, 
                                                             family = family_type, nfolds = ifelse(is.null(minnum) == FALSE, ceiling(nrow(cval) / minnum), 
                                                                                                   nfolds), type.measure = measure_type, alpha = alpha_val,   
                                                             gamma = seq(0, 1, length.out = ngamma), relax = relax_me, parallel = FALSE))
                            
                            
                            nntf <- nnt_find(fitb = cv.fitB, xmat = mf, vecy = yss, var = nnt_var, itm = nnt_tm, wtt = cval$missing_wt * wt)
                            m1b <- nntf$m1; m2b <- nntf$m2; m3b <- nntf$m3; m4b <- nntf$m4; m5b <- nntf$m5; 
                            p1b <- nntf$p1; p2b <- nntf$p2; p3b <- nntf$p3; p4b <- nntf$p4; p5b <- nntf$p5; 
                            inxb <- nntf$inx; tmb <- nntf$tm
                            
                            
                          } else if(psm == FALSE & nnt == TRUE){  
                            system.time(cv.fitB <- cv.glmnet(mf, yss, family = family_type, weights = cval$missing_wt * wt, 
                                                             nfolds = ifelse(is.null(minnum) == FALSE, ceiling(nrow(cval) / minnum), nfolds), type.measure = measure_type, alpha = alpha_val,   
                                                             gamma = seq(0, 1, length.out = ngamma), relax = relax_me, parallel = FALSE))
                            
                            nntf <- nnt_find(fitb = cv.fitB, xmat = mf, vecy = yss, var = nnt_var, itm = nnt_tm, wtt = cval$missing_wt * wt)
                            m1b <- nntf$m1; m2b <- nntf$m2; m3b <- nntf$m3; m4b <- nntf$m4; m5b <- nntf$m5; 
                            p1b <- nntf$p1; p2b <- nntf$p2; p3b <- nntf$p3; p4b <- nntf$p4; p5b <- nntf$p5; 
                            inxb <- nntf$inx; tmb <- nntf$tm
                            
                          }  else if(nnt == FALSE){
                            system.time(cv.fitB <- cv.glmnet(mf, yss, family = family_type, weights = cval$missing_wt * wt, 
                                                             nfolds = ifelse(is.null(minnum) == FALSE, ceiling(nrow(cval) / minnum), nfolds), type.measure = measure_type, alpha = alpha_val,   
                                                             gamma = seq(0, 1, length.out = ngamma), relax = relax_me, parallel = FALSE))
                            m1b <- m2b <- m3b <- m4b <- m5b <- NULL
                            
                          }
                          
                          a <- coef(cv.fitB, s = "lambda.min")
                          b <- m1b
                          c <- m2b
                          d <- m3b
                          e <- m4b
                          f <- m5b
                          
                          list(a, b, c, d, e, f)
                          
                        })
      
      stopCluster(cl)
      
      betaB <- rall[[1]]
      betaB <- do.call(cbind, lapply(seq_len(length(betaB)), function(bb) betaB[[bb]]))
      m1b <- rall[[2]]    
      m2b <- rall[[3]]     
      m3b <- rall[[4]]
      m4b <- rall[[5]]
      m5b <- rall[[6]]
      
      rn <- name_coef
      betaB <- betaB[name_coef2 %!in% grep("spline", name_coef2, value = T), ]
      
      ## bias-adjusted based on fractional-random weight bootstrap 
      boot_bound <- do.call(rbind, lapply(seq_len(length(beta.fit)), function(jj) {
        ic <- which(rn == rn[jj]) 
        icb <- which(rownames(betaB) == rn[jj])
        bb <- sort(betaB[icb, ])
        q <- length(which(bb < beta.fit[ic])) / B
        lb_ix <- round(B * pnorm((2 * qnorm(q, 0, 1)) + qnorm(alpha / 2, 0, 1)))
        lb_ix <- ifelse(lb_ix == 0, 1, lb_ix)
        ub_ix <- round(B * pnorm((2 * qnorm(q, 0, 1)) + qnorm(1 - (alpha / 2), 0, 1)))
        data.frame("names" =  rn[jj], "lb_bias2" = bb[lb_ix], "ub_bias2" = bb[ub_ix])
      }))
      
      
      
      bias <- (2 * as.numeric(beta.fit)) - apply(betaB, 1, mean) 
      dat <- data.frame("Variable" = name_coef,
                        "bOR" = exp(as.numeric(beta.fit)),
                        "bOR_050" = exp(apply(betaB, 1, finq050)),
                        "lb_025" = exp(apply(betaB, 1, finq025)),
                        "ub_975" = exp(apply(betaB, 1, finq975)),
                        "sd" = sqrt(apply(betaB, 1, var)),
                        
                        "p_val" =  2 * (1 - pnorm(abs(as.numeric(beta.fit) / sqrt(apply(betaB, 1, var))))), 
                        
                        "lb_bias" = exp(bias - ME * sqrt(apply(betaB, 1, var)) ),
                        "ub_bias" = exp(bias + ME * sqrt(apply(betaB, 1, var)) ),
                        
                        "lb_bias2" = exp(boot_bound$lb_bias2), 
                        "ub_bias2" = exp(boot_bound$ub_bias2),
                        
                        "lb_emp" = exp(2 * as.numeric(beta.fit) - apply(betaB, 1, finq975)),
                        "ub_emp" = exp(2 * as.numeric(beta.fit) + apply(betaB, 1, finq025)),
                        "lb" = exp(as.numeric(beta.fit) - ME * sqrt(apply(betaB, 1, var)) ),  
                        "ub" = exp(as.numeric(beta.fit) + ME * sqrt(apply(betaB, 1, var)) ))
      
    } 
  }
  
  # survival with imputation
  if(logistic == FALSE & imputeme == TRUE){
    
    dat_imp <- bm <- beta.fit <- trt.fitb <- rmst.fitb <- varimpb <- trt.fit <- rmst.fit <- varimp  <- list(); 
    
    tm <- inx <- m1 <- m2 <- m3 <- m4 <- m5 <- m1b <- m2b <- m3b <- m4b <- m5b <- p1 <- p2 <- p3 <- p4 <- p5 <- list()
    
    for(ll in 1 : nsample){ # run for each imputed dataset 
      
      var_to_add_adj <- nnt_var <- var_to_add[var_to_add %!in% mc$dropme]
      cval_imp <- mc$subd_imp[[ll]][, colnames(mc$subd_imp[[ll]]) %!in% mc$dropme]
      
      cval_imp$time <- cval_imp$time # scaling time
      
      nasum <- function(xx) {sum(is.na(xx))}
      
      if(length(which(apply(cval_imp[, colnames(cval_imp) %!in% c(grep("FUP30D_", colnames(cval_imp), value = T), 
                                                                  grep("OPTIMUM_", colnames(cval_imp), value = T),
                                                                  grep("_COMBO", colnames(cval), value = T),
                                                                  grep("ate_wt", colnames(cval_imp), value = T),
                                                                  grep("l_", colnames(cval_imp), value = T))], 2, nasum) >= 1)) >= 1){
        stop("Missing labs. Add unknown label")
      }
      
      if(psm == TRUE & psm_method == "splines"){
        library(splines2)
        
        try(l_vc <- cval_imp$l_vc)
        if(sum(is.na(l_vc)) != length(l_vc)){
          if(length(which(sort(l_vc) == max(sort(l_vc)))) > 1 | length(which(sort(l_vc) == min(sort(l_vc)))) > 1){ # adustment to by-pass error
            try(knots1 <- seq(head(sort(l_vc), 2)[2], tail(sort(l_vc), 2)[1] - 0.001, length.out = 6))  
          } else{
            try(knots1 <- seq(head(sort(l_vc), 2)[2], tail(sort(l_vc), 2)[1], length.out = 6))
          }
          if(min(knots1) <= min(l_vc)){
            knots1[1] <- knots1[1] + 0.001
          }
          if(max(knots1) <= max(l_vc)){
            knots1[length(knots1)] <- knots1[length(knots1)] - 0.001
          }
          try(bs_l_vc <- bSpline(l_vc, knots = knots1, degree = 0, intercept = TRUE))
        }
        
        try(l_vac <- cval_imp$l_vac)
        if(sum(is.na(l_vac)) != length(l_vac)){
          if(length(which(sort(l_vac) == max(sort(l_vac)))) > 1 | length(which(sort(l_vac) == min(sort(l_vac)))) > 1){ # adustment to by-pass error
            try(knots2 <- seq(head(sort(l_vac), 2)[2], tail(sort(l_vac), 2)[1] - 0.001, length.out = 6))  
          } else{
            try(knots2 <- seq(head(sort(l_vac), 2)[2], tail(sort(l_vac), 2)[1], length.out = 6))
          }
          if(min(knots2) <= min(l_vac)){
            knots2[1] <- knots2[1] + 0.001
          }
          if(max(knots2) <= max(l_vac)){
            knots2[length(knots2)] <- knots2[length(knots2)] - 0.001
          }
          
          try(bs_l_vac <- bSpline(l_vac, knots = knots2, degree = 0, intercept = TRUE))
        }
        
        try(l_a75c10 <- cval_imp$l_a75c10)
        if(sum(is.na(l_a75c10)) != length(l_a75c10)){
          if(length(which(sort(l_a75c10) == max(sort(l_a75c10)))) > 1 | length(which(sort(l_a75c10) == min(sort(l_a75c10)))) > 1){ # adustment to by-pass error
            try(knots3 <- seq(head(sort(l_a75c10), 2)[2], tail(sort(l_a75c10), 2)[1] - 0.001, length.out = 6))  
          } else{
            try(knots3 <- seq(head(sort(l_a75c10), 2)[2], tail(sort(l_a75c10), 2)[1], length.out = 6))
          }
          if(min(knots3) <= min(l_a75c10)){
            knots3[1] <- knots3[1] + 0.001
          }
          if(max(knots3) <= max(l_a75c10)){
            knots3[length(knots3)] <- knots3[length(knots3)] - 0.001
          }
          try(bs_l_a75c10 <- bSpline(l_a75c10, knots = knots3, degree = 0, intercept = TRUE))
        }
        
        try(l_a75c20 <- cval_imp$l_a75c20)
        if(sum(is.na(l_a75c20)) != length(l_a75c20)){
          if(length(which(sort(l_a75c20) == max(sort(l_a75c20)))) > 1 | length(which(sort(l_a75c20) == min(sort(l_a75c20)))) > 1){ # adustment to by-pass error
            try(knots4 <- seq(head(sort(l_a75c20), 2)[2], tail(sort(l_a75c20), 2)[1] - 0.001, length.out = 6))  
          } else{
            try(knots4 <- seq(head(sort(l_a75c20), 2)[2], tail(sort(l_a75c20), 2)[1], length.out = 6))
          }
          if(min(knots4) <= min(l_a75c20)){
            knots4[1] <- knots4[1] + 0.001
          }
          if(max(knots4) <= max(l_a75c20)){
            knots4[length(knots4)] <- knots4[length(knots4)] - 0.001
          }
          try(bs_l_a75c20 <- bSpline(l_a75c20, knots = knots4, degree = 0, intercept = TRUE))
        }
        
        
        if(alltx == "all") {
          bs_mat <- cbind(bs_l_vc, bs_l_vac, bs_l_a75c10, bs_l_a75c20)
        } else if(alltx == "a75c20 vs va"){ # drop vc
          bs_mat <- cbind(bs_l_a75c20)
        } else if(alltx == "a75c10 vs va"){ # drop vc
          bs_mat <- cbind(bs_l_a75c10)
        }
        
        colnames(bs_mat) <- paste0("spline", 1 : ncol(bs_mat), "")
      }
      
      nm0 <- paste0("", var_to_add_adj, "", collapse = " + " )
      part1 <- paste0("Surv(time, censor) ~ ", nm0, "", collapse = " + ")
      
      if(method == "Cox-PH + Penalty + Boot"){
        
        fom <- paste0("survival::coxph(", part1, ", ties = 'breslow', iter.max = 350, 
                  outer.max = 500, robust = TRUE,
                  singular.ok = TRUE, eps = 1e-1,
                  toler.chol = .Machine$double.eps^.35,
                  data = cval_imp)")
        try(fom <- eval(parse(text =  fom)))
        
        # extract model matrix
        mf <- mf_back <- model.matrix(fom)
        
        yss <- Surv(cval_imp$time, cval_imp$censor)
        
        set.seed(AA + 111 * ll)
        if(psm == TRUE & psm_method == "IPTW"){
          system.time(cv.fit <- cv.glmnet(mf, yss, weights = cval_imp$missing_wt * cval_imp$ate_wt1, 
                                          family = family_type, nfolds = ifelse(is.null(minnum) == FALSE, ceiling(nrow(cval_imp) / minnum), 
                                                                                nfolds), type.measure = measure_type, alpha = alpha_val,   
                                          gamma = seq(0, 1, length.out = ngamma), relax = relax_me, parallel = FALSE))
        } else if(psm == TRUE & psm_method == "splines"){
          mf <- cbind(mf_back, bs_mat)
          system.time(cv.fit <- cv.glmnet(mf, yss, weights = cval_imp$missing_wt, 
                                          family = family_type, nfolds = ifelse(is.null(minnum) == FALSE, ceiling(nrow(cval_imp) / minnum), 
                                                                                nfolds), type.measure = measure_type, alpha = alpha_val,   
                                          gamma = seq(0, 1, length.out = ngamma), relax = relax_me, parallel = FALSE))
        } else if(psm == FALSE){
          system.time(cv.fit <- cv.glmnet(mf, yss, weights = cval_imp$missing_wt, 
                                          family = family_type, nfolds = ifelse(is.null(minnum) == FALSE, ceiling(nrow(cval_imp) / minnum), 
                                                                                nfolds), type.measure = measure_type, alpha = alpha_val,   
                                          gamma = seq(0, 1, length.out = ngamma), relax = relax_me, parallel = FALSE))
        }  
        
        beta.fit[[ll]] <- coef(cv.fit, s = "lambda.min") # This is equivalent ot glmnet(mf, yss, family = "cox", lambda = cv.fit$lambda.min)
        name_coef <- name_coef2 <- rownames(beta.fit[[ll]])
        name_coef <- name_coef[name_coef2 %!in% grep("spline", name_coef2, value = T)]
        beta.fit[[ll]] <- beta.fit[[ll]][name_coef2 %!in% grep("spline", name_coef2, value = T), 1]
        time <- cval_imp$time
        censor <- cval_imp$censor 
        t.unique <- sort(unique(time[cval$censor == 1L]))
        centered = FALSE 
        
        
        if(nnt == TRUE){ 
          nntf <- nnt_find(fitb = cv.fit, xmat = mf, vecy = yss, var = nnt_var, itm = nnt_tm, wtt = cval$missing_wt)
          m1[[ll]] <- nntf$m1; m2[[ll]] <- nntf$m2; m3[[ll]] <- nntf$m3; m4[[ll]] <- nntf$m4; m5[[ll]] <- nntf$m5; 
          p1[[ll]] <- nntf$p1; p2[[ll]] <- nntf$p2; p3[[ll]] <- nntf$p3; p4[[ll]] <- nntf$p4; p5[[ll]] <- nntf$p5; 
          inx[[ll]] <- nntf$inx; 
          tm[[ll]] <- nntf$tm
        }
        
        # parallel computation 
        comb <- function(x, ...) {
          lapply(seq_along(x),
                 function(i) c(x[[i]], lapply(list(...), function(y) y[[i]])))
        }
        
        
        cores <- detectCores()
        cl <- makeCluster(cores[1] - 11) 
        registerDoParallel(cl)
        
        
        system.time(
          rall <- foreach(bb = 1 : B, .combine = 'comb', .multicombine=TRUE, .errorhandling = 'pass',
                       
                          .init = list(list(), list(), list(), list(), list(), list()), 
                          .export = c("nnt_find",  '%!in%'), 
                          .packages = c('glmnet',  'ranger', 'doParallel', 'survival')) %dopar% {
                            
                            unregister <- function() {
                              env <- foreach:::.foreachGlobals
                              rm(list=ls(name=env), pos=env)
                            }    
                            
                            unregister()
                            print(paste0("now-", bb, ""))     
                            
                            set.seed(bb + AA + 111)
                            n <- nrow(cval_imp)
                            wt <- rexp(n, 1)
                            wt <- wt / sum(wt)
                          
                            
                            if(psm == TRUE & psm_method == "IPTW" & nnt == TRUE){
                              system.time(cv.fitB <- cv.glmnet(mf, yss, weights = cval_imp$missing_wt * wt * cval_imp$ate_wt1, 
                                                               family = family_type, nfolds = ifelse(is.null(minnum) == FALSE, ceiling(nrow(cval_imp) / minnum), 
                                                                                                     nfolds), type.measure = measure_type, alpha = alpha_val,   
                                                               gamma = seq(0, 1, length.out = ngamma), relax = relax_me, parallel = FALSE))
                              
                              nntf <- nnt_find(fitb = cv.fitB, xmat = mf, vecy = yss, var = nnt_var, itm = nnt_tm, wtt = cval$missing_wt * wt * cval$ate_wt1)
                              m1b2 <- nntf$m1; m2b2 <- nntf$m2; m3b2 <- nntf$m3; m4b2 <- nntf$m4; m5b2 <- nntf$m5; 
                              p1b <- nntf$p1; p2b <- nntf$p2; p3b <- nntf$p3; p4b <- nntf$p4; p5b <- nntf$p5; 
                              inxb <- nntf$inx; tmb <- nntf$tm
                              
                            } else if(psm == TRUE & psm_method == "splines" & nnt == TRUE){
                              mf <- cbind(mf_back, bs_mat)
                              system.time(cv.fitB <- cv.glmnet(mf, yss, weights = cval_imp$missing_wt * wt, 
                                                               family = family_type, nfolds = ifelse(is.null(minnum) == FALSE, ceiling(nrow(cval_imp) / minnum), 
                                                                                                     nfolds), type.measure = measure_type, alpha = alpha_val,   
                                                               gamma = seq(0, 1, length.out = ngamma), relax = relax_me, parallel = FALSE))
                              
                              if(is.null(cval$ate_wt1) == TRUE){
                                nntf <- nnt_find(fitb = cv.fitB, xmat = mf, vecy = yss, var = nnt_var, itm = nnt_tm, wtt = cval$missing_wt * wt)
                              } else{
                                nntf <- nnt_find(fitb = cv.fitB, xmat = mf, vecy = yss, var = nnt_var, itm = nnt_tm, wtt = cval$missing_wt * wt * cval$ate_wt1)
                              }
                              m1b2 <- nntf$m1; m2b2 <- nntf$m2; m3b2 <- nntf$m3; m4b2 <- nntf$m4; m5b2 <- nntf$m5; 
                              p1b <- nntf$p1; p2b <- nntf$p2; p3b <- nntf$p3; p4b <- nntf$p4; p5b <- nntf$p5; 
                              inxb <- nntf$inx; tmb <- nntf$tm
                              
                            } else if(psm == FALSE & nnt == TRUE){  
                              system.time(cv.fitB <- cv.glmnet(mf, yss, family = family_type, weights = cval_imp$missing_wt * wt, 
                                                               nfolds = ifelse(is.null(minnum) == FALSE, ceiling(nrow(cval_imp) / minnum), nfolds), type.measure = measure_type, alpha = alpha_val,   
                                                               gamma = seq(0, 1, length.out = ngamma), relax = relax_me, parallel = FALSE))
                              
                              if(is.null(cval$ate_wt1) == TRUE){
                                nntf <- nnt_find(fitb = cv.fitB, xmat = mf, vecy = yss, var = nnt_var, itm = nnt_tm, wtt = cval$missing_wt * wt)
                              } else{
                                nntf <- nnt_find(fitb = cv.fitB, xmat = mf, vecy = yss, var = nnt_var, itm = nnt_tm, wtt = cval$missing_wt * wt * cval$ate_wt1)
                              }
                              m1b2 <- nntf$m1; m2b2 <- nntf$m2; m3b2 <- nntf$m3; m4b2 <- nntf$m4; m5b2 <- nntf$m5; 
                              p1b <- nntf$p1; p2b <- nntf$p2; p3b <- nntf$p3; p4b <- nntf$p4; p5b <- nntf$p5; 
                              inxb <- nntf$inx; tmb <- nntf$tm
                            } else if(nnt == FALSE){
                              system.time(cv.fitB <- cv.glmnet(mf, yss, family = family_type, weights = cval$missing_wt * wt, 
                                                               nfolds = ifelse(is.null(minnum) == FALSE, ceiling(nrow(cval) / minnum), nfolds), type.measure = measure_type, alpha = alpha_val,   
                                                               gamma = seq(0, 1, length.out = ngamma), relax = relax_me, parallel = FALSE))
                              m1b2 <- m2b2 <- m3b2 <- m4b2 <- m5b2 <- NULL
                              
                            }
                            
                            a <- coef(cv.fitB, s = "lambda.min")
                            b <- m1b2
                            c <- m2b2
                            d <- m3b2
                            e <- m4b2
                            f <- m5b2
                            
                            list(a, b, c, d, e, f)
                            
                          })
        
        stopCluster(cl)
        
        betaB <- rall[[1]]
        betaB <- do.call(cbind, lapply(seq_len(length(betaB)), function(bb) betaB[[bb]]))
        m1b[[ll]] <- rall[[2]]    
        m2b[[ll]] <- rall[[3]]     
        m3b[[ll]] <- rall[[4]]
        m4b[[ll]] <- rall[[5]]
        m5b[[ll]] <- rall[[6]]
        
        
        rn <- name_coef
        betaB <- betaB[name_coef2 %!in% grep("spline", name_coef2, value = T), ]
        
        ## bias adjusted based on fractional-random weight bootstrap https://par.nsf.gov/servlets/purl/10155761
        boot_bound <- do.call(rbind, lapply(seq_len(length(beta.fit[[ll]])), function(jj) {
          ic <- which(rn == rn[jj]) 
          icb <- which(rownames(betaB) == rn[jj])
          bb <- sort(betaB[icb, ])
          q <- length(which(bb < beta.fit[[ll]][ic])) / B
          lb_ix <- round(B * pnorm((2 * qnorm(q, 0, 1)) + qnorm(alpha / 2, 0, 1)))
          lb_ix <- ifelse(lb_ix == 0, 1, lb_ix)
          ub_ix <- round(B * pnorm((2 * qnorm(q, 0, 1)) + qnorm(1 - (alpha / 2), 0, 1)))
        
          data.frame("names" =  rn[jj], "lb_bias2" = bb[lb_ix], "ub_bias2" = bb[ub_ix])
        }))
        
        
        
        bias <- (2 * as.numeric(beta.fit[[ll]])) - apply(betaB, 1, mean) # following https://math.montana.edu/jobo/thainp/boot.pdf
      
        bm[[ll]] <- betaB
        
        
        dat_imp[[ll]] <- data.frame("Variable" = name_coef,
                                    "bOR" = exp(as.numeric(beta.fit[[ll]])),
                                    "bOR_050" = exp(apply(betaB, 1, finq050)),
                                    "lb_025" = exp(apply(betaB, 1, finq025)),
                                    "ub_975" = exp(apply(betaB, 1, finq975)),
                                    "sd" = sqrt(apply(betaB, 1, var)),
                                    "p_val" =  2 * (1 - pnorm(abs(as.numeric(beta.fit[[ll]]) / sqrt(apply(betaB, 1, var))))),
                                    "lb_bias" = exp(bias - ME * sqrt(apply(betaB, 1, var)) ),
                                    "ub_bias" = exp(bias + ME * sqrt(apply(betaB, 1, var)) ),
                                    
                                    "lb_bias2" = exp(boot_bound$lb_bias2), 
                                    "ub_bias2" = exp(boot_bound$ub_bias2),
                                    
                                    "lb_emp" = exp(2 * as.numeric(beta.fit[[ll]]) - apply(betaB, 1, finq975)),
                                    "ub_emp" = exp(2 * as.numeric(beta.fit[[ll]]) + apply(betaB, 1, finq025)),
                                   
                                    "lb" = exp(as.numeric(beta.fit[[ll]]) - ME * sqrt(apply(betaB, 1, var)) ),  
                                    "ub" = exp(as.numeric(beta.fit[[ll]]) + ME * sqrt(apply(betaB, 1, var)) ))
       
      } # methods ends 
      
      
      ######## Impute == TRUE 
      if(method %in% c("CausalML", "CausalML1000", "CausalML120", "CausalML365")){
        dat <- NULL 
        nm0 <- paste0("", var_to_add[var_to_add %!in% "TX_TYPE_CML"], "", collapse = " + " )
        part1 <- paste0("Surv(time, censor) ~ ", nm0, "", collapse = " + ")
        
        
        fom <- paste0("survival::coxph(", part1, ", ties = 'breslow', iter.max = 1000, 
                  outer.max = 500, robust = FALSE,
                  singular.ok = TRUE, eps = 1e-1,
                  toler.chol = .Machine$double.eps^.25,
                  data = cval_imp)")
        try(fom <- eval(parse(text =  fom)))
        
        # extract model matrix
        mf <- mf_back <- model.matrix(fom)
        dim(mf) 
        
        yss <- Surv(cval_imp$time, cval_imp$censor)
        X <- mf
        Y <- cval_imp$time 
        W <- ifelse(cval_imp$TX_TYPE_CML == "TRT", 1, 0)
        D <- cval_imp$censor
        table(W)
        require(grf)
        cs.forest11 <- causal_survival_forest(X, Y, W, D, horizon = c(30), target = "survival.probability")
        cs.forest12 <- causal_survival_forest(X, Y, W, D, horizon = c(90), target = "survival.probability")
        cs.forest13 <- causal_survival_forest(X, Y, W, D, horizon = c(180), target = "survival.probability")
        cs.forest14 <- causal_survival_forest(X, Y, W, D, horizon = c(365), target = "survival.probability")
        cs.forest15 <- causal_survival_forest(X, Y, W, D, horizon = c(450), target = "survival.probability")
        cs.forest16 <- causal_survival_forest(X, Y, W, D, horizon = c(730), target = "survival.probability")
        
        cs.forest21 <- causal_survival_forest(X, Y, W, D, horizon = c(30), target = "RMST")
        cs.forest22 <- causal_survival_forest(X, Y, W, D, horizon = c(90), target = "RMST")
        cs.forest23 <- causal_survival_forest(X, Y, W, D, horizon = c(180), target = "RMST")
        cs.forest24 <- causal_survival_forest(X, Y, W, D, horizon = c(365), target = "RMST")
        cs.forest25 <- causal_survival_forest(X, Y, W, D, horizon = c(450), target = "RMST")
        cs.forest26 <- causal_survival_forest(X, Y, W, D, horizon = c(730), target = "RMST")
        
        
        trt.fit[[ll]] <- rbind(average_treatment_effect(cs.forest11), average_treatment_effect(cs.forest12), 
                               average_treatment_effect(cs.forest13), average_treatment_effect(cs.forest14), 
                               average_treatment_effect(cs.forest15), average_treatment_effect(cs.forest16))
        rmst.fit[[ll]] <- rbind(average_treatment_effect(cs.forest21), average_treatment_effect(cs.forest22), 
                                average_treatment_effect(cs.forest23), average_treatment_effect(cs.forest24), 
                                average_treatment_effect(cs.forest25), average_treatment_effect(cs.forest26))
        
        #best_linear_projection(cs.forest, X)
        
        varimp[[ll]] <- cbind(variable_importance(cs.forest11), variable_importance(cs.forest12), 
                              variable_importance(cs.forest13), variable_importance(cs.forest14), 
                              variable_importance(cs.forest15), variable_importance(cs.forest16))
        rownames(varimp) <- colnames(X)
        
        # run below for parallal computation 
        comb <- function(x, ...) {
          lapply(seq_along(x),
                 function(i) c(x[[i]], lapply(list(...), function(y) y[[i]])))
        }
        
        cores <- detectCores()
        cl <- makeCluster(cores[1] - 2) #not to overload your computer
        registerDoParallel(cl)
        
        
        system.time(
          betaB <- foreach(i = seq_along(1 : B), .combine = 'comb', .multicombine=TRUE, .errorhandling = 'pass',
                           .init=list(list(), list(), list(), list()), 
                           .packages = c('glmnet',  'ranger', 'doParallel', 'survival', 'grf')) %dopar% {
                             #try(for(bb in 1 : B){
                             set.seed(i + 78787890 + AA)
                             n <- nrow(cval_imp)
                             wt <- rexp(n, 1)
                             wt <- wt / sum(wt)
                             
                             cs.forest11b <- causal_survival_forest(X, Y, W, D, horizon = c(30), target = "survival.probability", sample.weights = wt)
                             cs.forest12b <- causal_survival_forest(X, Y, W, D, horizon = c(90), target = "survival.probability", sample.weights = wt)
                             cs.forest13b <- causal_survival_forest(X, Y, W, D, horizon = c(180), target = "survival.probability", sample.weights = wt)
                             cs.forest14b <- causal_survival_forest(X, Y, W, D, horizon = c(365), target = "survival.probability", sample.weights = wt)
                             cs.forest15b <- causal_survival_forest(X, Y, W, D, horizon = c(450), target = "survival.probability", sample.weights = wt)
                             cs.forest16b <- causal_survival_forest(X, Y, W, D, horizon = c(730), target = "survival.probability", sample.weights = wt)
                             
                             cs.forest21b <- causal_survival_forest(X, Y, W, D, horizon = c(30), target = "RMST", sample.weights = wt)
                             cs.forest22b <- causal_survival_forest(X, Y, W, D, horizon = c(90), target = "RMST", sample.weights = wt)
                             cs.forest23b <- causal_survival_forest(X, Y, W, D, horizon = c(180), target = "RMST", sample.weights = wt)
                             cs.forest24b <- causal_survival_forest(X, Y, W, D, horizon = c(365), target = "RMST", sample.weights = wt)
                             cs.forest25b <- causal_survival_forest(X, Y, W, D, horizon = c(450), target = "RMST", sample.weights = wt)
                             cs.forest26b <- causal_survival_forest(X, Y, W, D, horizon = c(730), target = "RMST", sample.weights = wt)
                             
                             trt.fitb <- rbind(average_treatment_effect(cs.forest11b), average_treatment_effect(cs.forest12b), 
                                               average_treatment_effect(cs.forest13b), average_treatment_effect(cs.forest14b), 
                                               average_treatment_effect(cs.forest15b), average_treatment_effect(cs.forest16b))
                             rmst.fitb <- rbind(average_treatment_effect(cs.forest21b), average_treatment_effect(cs.forest22b), 
                                                average_treatment_effect(cs.forest23b), average_treatment_effect(cs.forest24b), 
                                                average_treatment_effect(cs.forest25b), average_treatment_effect(cs.forest26b))
                             
                             #best_linear_projection(cs.forest, X)
                             
                             varimpb <- cbind(variable_importance(cs.forest11b), variable_importance(cs.forest12b), 
                                              variable_importance(cs.forest13b), variable_importance(cs.forest14b), 
                                              variable_importance(cs.forest15b), variable_importance(cs.forest16b))
                             rownames(varimpb) <- colnames(X)
                             
                             A1 <- trt.fitb
                             B1 <- rmst.fitb
                             C1 <- varimpb
                             
                             list(A1, B1, C1, i)
                           })
        
        stopCluster(cl)
        
        
        trt.fitb[[ll]] <- betaB[[1]]
        rmst.fitb[[ll]] <- betaB[[2]]     # predicted probabilities of all (NOTE: It gives preidcted prob for both test and train data and first row is test data)
        varimpb[[ll]] <- betaB[[3]]     #  time of interest
        
      }
      ####### Impute === TRUE ends 
      
      
      
      
      
    }  # for each imp dataset 
    
    var_to_add <- var_to_add_adj
    
    
    ## biase based on fractional-random weight bootstrap https://par.nsf.gov/servlets/purl/10155761
    if(method %!in% c("CausalML", "CausalML1000", "CausalML120", "CausalML365")){
      rn <- rownames(beta.fit[[ll]])
      
      gt <- paste0("beta.fit[[", 1 : nsample, "]]", collapse  = " + ")
      gt <- eval(parse(text = gt))
      beta.fit2 <- gt / nsample
      
      dat <- data.frame("Variable" = name_coef,
                        
                        "bOR" = exp(as.numeric(beta.fit2)),
                        
                        
                        "bOR_050" = exp(unlist(lapply(seq_len(length(name_coef)), function(jj) { 
                          finq050(unlist(lapply(seq_len(nsample), function(ll) {  bm[[ll]][jj, ] })))
                        }))),
                        
                        
                        "lb_025" = exp(unlist(lapply(seq_len(length(name_coef)), function(jj) { 
                          finq025(unlist(lapply(seq_len(nsample), function(ll) {  bm[[ll]][jj, ] })))
                        }))),
                        
                        "ub_975" = exp(unlist(lapply(seq_len(length(name_coef)), function(jj) { 
                          finq975(unlist(lapply(seq_len(nsample), function(ll) {  bm[[ll]][jj, ] })))
                        }))))
      dat$lb_bias2 <- dat$lb <- dat$lb_025 
      dat$ub_bias2 <- dat$ub <- dat$ub_975 
      
      if(nnt == TRUE & imputeme == TRUE){
        m1r <- m2r <- m3r <- m4r <- m5r <- g1b <- g2b <- g3b <- g4b <- g5b <- list()
        
        for(rr in 1 : length(m1[[1]])){
          
          gt <- paste0("m1[[", 1 : nsample, "]][[rr]]", collapse  = " + ")
          gt <- eval(parse(text = gt))
          m1r[[rr]] <- gt / nsample
          
          
          
          gt <- paste0("m2[[", 1 : nsample, "]][[rr]]", collapse  = " + ")
          gt <- eval(parse(text = gt))
          m2r[[rr]] <- gt / nsample
          
          
          
          
          if(length(m3[[1]]) > 0 & imputeme == FALSE){
            gt <- paste0("m3[[", 1 : nsample, "]][[rr]]", collapse  = " + ")
            gt <- eval(parse(text = gt))
            m3r[[rr]] <- gt / nsample
            
            
            
          }
          
          if(length(m4[[1]]) > 0 & imputeme == FALSE){
            gt <- paste0("m4[[", 1 : nsample, "]][[rr]]", collapse  = " + ")
            gt <- eval(parse(text = gt))
            m4r[[rr]] <- gt / nsample
            
            
          }
          
          if(length(m5[[1]]) > 0 & imputeme == FALSE){
            gt <- paste0("m5[[", 1 : nsample, "]][[rr]]", collapse  = " + ")
            gt <- eval(parse(text = gt))
            m5r[[rr]] <- gt / nsample
            
            
            
            
          }
          
          
        }
        
        m1 <- m1r;  m2 <- m2r;  m3 <- m3r;  m4 <- m4r;  m5 <- m5r
        
        # collapsing over nsamples for bootstrap values 
        QQ <- paste0("try(g1b <- lapply(seq_len(B), function(uu) {
          kk <- lapply(seq_len(length(m1)), function(vv) {
           (m1b[[1]][[uu]][[vv]] + m1b[[2]][[uu]][[vv]] + m1b[[3]][[uu]][[vv]] + m1b[[4]][[uu]][[vv]] + m1b[[5]][[uu]][[vv]]) * (1 / 5)
           })
           kk
          }))") 
        QQ <- eval(parse(text = QQ))
        
        # collapsing over nsamples for bootstrap values 
        QQ <- paste0("try(g2b <- lapply(seq_len(B), function(uu) {
          kk <- lapply(seq_len(length(m2)), function(vv) {
           (m2b[[1]][[uu]][[vv]] + m2b[[2]][[uu]][[vv]] + m2b[[3]][[uu]][[vv]] + m2b[[4]][[uu]][[vv]] + m2b[[5]][[uu]][[vv]]) * (1 / 5)
           })
           kk
          }))") 
        QQ <- eval(parse(text = QQ))
        
        # collapsing over nsamples for bootstrap values 
        if(length(m3) > 0 & imputeme == FALSE){
          QQ <- paste0("try(g3b <- lapply(seq_len(B), function(uu) {
          kk <- lapply(seq_len(length(m3)), function(vv) {
           (m3b[[1]][[uu]][[vv]] + m3b[[2]][[uu]][[vv]] + m3b[[3]][[uu]][[vv]] + m3b[[4]][[uu]][[vv]] + m3b[[5]][[uu]][[vv]]) * (1 / 5)
           })
           kk
          }))") 
          QQ <- eval(parse(text = QQ))
        }
        
        # collapsing over nsamples for bootstrap values 
        if(length(m4) > 0 & imputeme == FALSE){
          QQ <- paste0("try(g4b <- lapply(seq_len(B), function(uu) {
          kk <- lapply(seq_len(length(m4)), function(vv) {
           (m4b[[1]][[uu]][[vv]] + m4b[[2]][[uu]][[vv]] + m4b[[3]][[uu]][[vv]] + m4b[[4]][[uu]][[vv]] + m4b[[5]][[uu]][[vv]]) * (1 / 5)
           })
           kk
          }))") 
          QQ <- eval(parse(text = QQ))
        }
        
        
        # collapsing over nsamples for bootstrap values 
        if(length(m5) > 0 & imputeme == FALSE){
          QQ <- paste0("try(g5b <- lapply(seq_len(B), function(uu) {
          kk <- lapply(seq_len(length(m5)), function(vv) {
           (m5b[[1]][[uu]][[vv]] + m5b[[2]][[uu]][[vv]] + m5b[[3]][[uu]][[vv]] + m5b[[4]][[uu]][[vv]] + m5b[[5]][[uu]][[vv]]) * (1 / 5)
           })
           kk
          }))") 
          QQ <- eval(parse(text = QQ))
        }
        
        m1b <- g1b; m2b <- g2b; m3b <- g3b; m4b <- g4b; m5b <- g5b
        
        
      }
    }
    
    if(method %in% c("CausalML", "CausalML1000", "CausalML120", "CausalML365")){
      dat <- list(trt.fit = trt.fit, 
                  trt.fitb = trt.fitb, 
                  rmst.fit = rmst.fit, 
                  rmst.fitb = rmst.fitb, 
                  varimp = varimp, 
                  varimpb = varimpb)
    }
    
  } # loop closes
  
  dl <- NULL
  
  if(imputeme ==  FALSE){
    cval <- cval %>% mutate_if(is.factor, as.character)
  }
  
  if(imputeme == TRUE){
    cval <- cval_imp
    cval <- cval %>% mutate_if(is.character, as.factor)
  }
  
  
  if(delete_NP_param == TRUE){
    ## Start here deleting NP
    dl <- unlist(lapply(seq_len(nrow(dat)), function(ii) {
      length(grep("Missing", as.character(dat$Variable[ii]), value = F))}))
    ixx <- which(dl >= 1)
    if(length(ixx) > 0){
      dat <- dat[-ixx, ]
    }
    dat <- dat %>% mutate_if(is.factor, as.character)
    dat <- dat %>% mutate_if(is.character, as.factor)
  }
  
  
  output = list (dat = dat, var_to_add = var_to_add, data = cval,
                 dl_index = dl, 
                 tm = tm, inx = inx, 
                 m1 = m1, m2 = m2, m3 = m3, m4 = m4, m5 = m5,
                 m1b = m1b, m2b = m2b, m3b = m3b, m4b = m4b, m5b = m5b)
}





## Generate risk stratification 
# mcatg = A list containing covariate-level classification under different settings 
# fc_include = Boolean operator whether phenotypic features need to be included in risk model 
# os_type = Whether OS for both risk model and ELN is needed with pairwise comparison. ""os (ELN)" will not run pairwise comparison
# mis_type = Covariate-level classification features were estimated based on what missing data model (IPW, Category, MICE)
# cen_type = indicator how allo-HCT patients were treated
# global_var = Non-null value indicates practical model highlighting only features that are based on CYT and composite mutation (used for the RRM model generation; https://www.medrxiv.org/content/10.1101/2024.12.02.24318344v1)
# good_def = Whether neutral positive variables to be added in Favorable risk group 
# exclude_cd71 = "No" is associated with Rule-I and Rule-III and "5PY" with Rule-II and Rule-IV
# data_type = "predictive (cv)" corresponds to evaluating predictive performances of risk model. "non-predictive" does not run this. 
# test_ix = IDs for test subject 
# t_max = Maximum time in the dynamic risk profile 
# tux = Time-points at which marginal risks need to be evaluated
# vrm = Covariates to be adjusted in predictive model 
# getd = Directory link where analytical dataset is saved. If there is a test dataset, both training and test datasets should be concatenated row-wise.  
findkm_riskgroup_paperML <- function(mcatg, fc_include = TRUE, os_type = "os (ELN) + pairwise", mis_type = "catg", cen_type = "all", 
                                     global_var = NULL, good_def = "yes_neutral", exclude_cd71 = "No", data_type = "non-predictive", 
                                     test_ix = NULL, t_max = 365 * 4, tux = seq(30, 450, by = 30), 
                                     vrm = c("DEM_AGE_RANGE", "DEM_RACE", "DEM_SEX", "MOD_RISK", "MOD_RISK_EDIT"), 
                                             getd = "/mnt/filestore/bkt-dev-prj-dev-rfsci-shiny-nonprt-naz/illustration/rds"){
 
  QQ <- paste0("setwd('", getd, "/')"); QQ <- eval(parse(text = QQ))
  QQ <- paste0("KMdata <- KMdata_va <- readRDS('", getd, "/KMdata.rds')"); QQ <- eval(parse(text = QQ))
  QQ <- paste0("KMdata_sct <- KMdata_va_sct <- readRDS('", getd, "/KMdata_sct.rds')"); QQ <- eval(parse(text = QQ))
    
  
  if(fc_include == TRUE){  
    mcatg <- mcatg
    
  } else if(fc_include == FALSE){
    mcatg <- mcatg[-grep("FC_", mcatg$Variable, value = F), ]
    
  }
  
  sa <- as.character(mcatg$Variable[which(mcatg$Values == "Strongly adverse")])
  ma <- as.character(mcatg$Variable[which(mcatg$Values == "Moderately adverse")])
  na <- as.character(mcatg$Variable[which(mcatg$Values == "Neutral signaling adversity")])
  mf <- as.character(mcatg$Variable[which(mcatg$Values == "Moderately favorable")])
  sf <- as.character(mcatg$Variable[which(mcatg$Values == "Strongly favorable")])
  nf <- as.character(mcatg$Variable[which(mcatg$Values == "Neutral signaling favorability")])
  
  sa_edit <- as.character(mcatg$Variable[which(mcatg$Values_mean == "Strongly adverse")])
  ma_edit <- as.character(mcatg$Variable[which(mcatg$Values_mean == "Moderately adverse")])
  na_edit <- as.character(mcatg$Variable[which(mcatg$Values_mean == "Neutral signaling adversity")])
  mf_edit <- as.character(mcatg$Variable[which(mcatg$Values_mean == "Moderately favorable")])
  sf_edit <- as.character(mcatg$Variable[which(mcatg$Values_mean == "Strongly favorable")])
  nf_edit <- as.character(mcatg$Variable[which(mcatg$Values_mean == "Neutral signaling favorability")])
  
  
  if(exclude_cd71 == "5PY"){ # excluding variables that were too few or too prevalent in CU
    
    ltr <- c(sa, ma, sf, mf, na, nf)
    ltr_ex <- unlist(lapply(seq_len(length(ltr)), function(yy) {
      
      FF1 <- paste0("sumy(KMdata_va$",  ltr[yy], ")")
      FF1 <- eval(parse(text = FF1))
      
      FF2 <- paste0("sumn(KMdata_va$",  ltr[yy], ")")
      FF2 <- eval(parse(text = FF2))
      
      FF <- FF1 / (FF1 + FF2)
      
    }))
    
    vl <- 0.40
    # deleting variables that are too prevalent, variables with selection bias, or too few cases
    sa <- sa[sa %!in% ltr[which(ltr_ex < 0.02 | ltr_ex > vl)]]
    ma <- ma[ma %!in% ltr[which(ltr_ex < 0.02 | ltr_ex > vl)]]
    
    sa_edit <- sa_edit[sa_edit %!in% ltr[which(ltr_ex < 0.02 | ltr_ex > vl)]]
    ma_edit <- ma_edit[ma_edit %!in% ltr[which(ltr_ex < 0.02 | ltr_ex > vl)]]
    
    sf <- sf[sf %!in% ltr[which(ltr_ex < 0.02 | ltr_ex > vl)]]
    mf <- mf[mf %!in% ltr[which(ltr_ex < 0.02 | ltr_ex > vl)]]
    
    sf_edit <- sf_edit[sf_edit %!in% ltr[which(ltr_ex < 0.02 | ltr_ex > vl)]]
    mf_edit <- mf_edit[mf_edit %!in% ltr[which(ltr_ex < 0.02 | ltr_ex > vl)]]
    
  }
  
  sa_back <- sa; sf_back <- sf; mf_back <- mf; ma_back <- ma; na_back <- na; nf_back <- nf
  sa_edit_back <- sa_edit; sf_edit_back <- sf_edit; mf_edit_back <- mf_edit; ma_edit_back <- ma_edit; na_edit_back <- na_edit; nf_edit_back <- nf_edit
  
  yesme <- function(gg) {length(which(gg == "Yes"))}
  
  
  # original data
  KMdatas <- KMdata[, colnames(KMdata) %in% c(sa, ma, mf, sf, na, nf, sa_edit, ma_edit, mf_edit, sf_edit, na_edit, nf_edit, "time", "censor",   grep("_MUTATION", colnames(KMdata), value = T))] 
  
  
  
  
  ## checking where a subject belongs to
  wght <- wght_fin <- do.call(rbind, lapply(seq_len(nrow(KMdatas)), function(ii) {
    
    data.frame("bad" = yesme(KMdatas[ii, colnames(KMdatas) %in% c(sa, ma)]), 
               "neutral_bad" = yesme(KMdatas[ii, colnames(KMdatas) %in% c(na)]), 
               "good" = yesme(KMdatas[ii, colnames(KMdatas) %in% c(sf, mf)]), 
               "neutral_good" = yesme(KMdatas[ii, colnames(KMdatas) %in% c(nf)]))
  }))
  
  wght_edit <- wght_fin_edit <- do.call(rbind, lapply(seq_len(nrow(KMdatas)), function(ii) {
    
    data.frame("bad" = yesme(KMdatas[ii, colnames(KMdatas) %in% c(sa_edit, ma_edit)]), 
               "neutral_bad" = yesme(KMdatas[ii, colnames(KMdatas) %in% c(na_edit)]), 
               "good" = yesme(KMdatas[ii, colnames(KMdatas) %in% c(sf_edit, mf_edit)]), 
               "neutral_good" = yesme(KMdatas[ii, colnames(KMdatas) %in% c(nf_edit)]))
  }))
  
  
  
  adv_ix <- which(wght$bad > 0 & wght$good == 0)
  int_ix <- which(wght$bad > 0 & wght$good > 0)
  int_fav_ix <- which(wght$bad == 0 & wght$good == 0 & wght$neutral_good > 0)
  fav_ix <- which(wght$bad == 0 & wght$good > 0 )
  neu_ix <- which(wght$bad == 0 & wght$good == 0 & wght$neutral_good == 0)
  
  
  adv_edit_ix <- which(wght_edit$bad > 0 & wght_edit$good == 0)
  int_edit_ix <- which(wght_edit$bad > 0 & wght_edit$good > 0)
  int_edit_fav_ix <- which(wght_edit$bad == 0 & wght_edit$good == 0 & wght_edit$neutral_good > 0)
  fav_edit_ix <- which(wght_edit$bad == 0 & wght_edit$good > 0 )
  neu_edit_ix <- which(wght_edit$bad == 0 & wght_edit$good == 0 & wght_edit$neutral_good == 0)
  
  
  ## 1st def
  KMdata$MOD_RISK <- rep("Intermediate", nrow(KMdata))
  KMdata$MOD_RISK[adv_ix] <- "Adverse"
  KMdata$MOD_RISK[fav_ix] <- "Favorable"
  table(KMdata$MOD_RISK)
  
  if(good_def == "yes_neutral"){
    KMdata$MOD_RISK[int_fav_ix] <- "Favorable"
  } 
  KMdata_test <- KMdata
  table(KMdata$MOD_RISK)
  
  if("Favorable" %in% names(table(KMdata$MOD_RISK))){
    try(KMdata <- KMdata_test <- KMdata %>% mutate(MOD_RISK = fct_relevel(MOD_RISK, c("Favorable", "Intermediate", "Adverse")), 
                                                   MOD_RISK = relevel(MOD_RISK, ref = 'Favorable')));
  }
  
  ## 2nd def
  KMdata$MOD_RISK_EDIT <- rep("Intermediate", nrow(KMdata))
  KMdata$MOD_RISK_EDIT[adv_edit_ix] <- "Adverse"
  KMdata$MOD_RISK_EDIT[fav_edit_ix] <- "Favorable"
  if(good_def == "yes_neutral"){
    KMdata$MOD_RISK_EDIT[int_fav_edit_ix] <- "Favorable"
  } 
  KMdata_test <- KMdata
  table(KMdata$MOD_RISK_EDIT)
  if("Favorable" %in% names(table(KMdata$MOD_RISK_EDIT))){
    try(KMdata <- KMdata_test <- KMdata %>% mutate(MOD_RISK_EDIT = fct_relevel(MOD_RISK_EDIT, c("Favorable", "Intermediate", "Adverse")), 
                                                   MOD_RISK_EDIT = relevel(MOD_RISK_EDIT, ref = 'Favorable')));
  }
  
  
  KMdata_sct$MOD_RISK <- KMdata$MOD_RISK
  KMdata_sct$MOD_RISK_EDIT <- KMdata$MOD_RISK_EDIT
  
  
  
  if(is.null(global_var) == FALSE) {
    
    ## 1st def
    sa <- sa[sa %in% c(grep("CYT_", sa, value = T), grep("TP53_", sa, value = T), grep("CYT_MINUS7", sa, value = T), grep("CYT_MINUS5", sa, value = T))]
    ma <- ma[ma %in% c(grep("CYT_", ma, value = T), grep("TP53_", ma, value = T), grep("CYT_MINUS7", ma, value = T), grep("CYT_MINUS5", ma, value = T))] # c("CYT_DEL7Q", "TP53_MUTATION", "CYT_COMPLEX", "CYT_INV3")
    
    sf <- sf[sf %in% c(grep("CYT_", sf, value = T), grep("IDH1", sf, value = T), grep("IDH2", sf, value = T), grep("NPM1", sf, value = T))] # c( "IDH1_MUTATION", "IDH2_MUTATION", "NPM1_MUTATION", "NGS_IDH1", "NGS_IDH2", "NGS_NPM1")
    mf <- mf[mf %in% c(grep("CYT_", mf, value = T), grep("IDH1", mf, value = T), grep("IDH2", mf, value = T), grep("NPM1", mf, value = T))] # c( "IDH1_MUTATION", "IDH2_MUTATION", "NPM1_MUTATION", "NGS_IDH1", "NGS_IDH2", "NGS_NPM1")
    
    na <- "FLT3ITD_MUTATION" #c("NGS_FLT3ITD", "PCR_FLT3_ITD")
    
    KMdatas <- KMdata[, colnames(KMdata) %in% c(sa, ma, mf, sf, na, nf, "time", "censor", "ELN_RISK_GROUP",
                                                "TP53_MUTATION", "NGS_KRAS", "NGS_NRAS")]
    
    dim(KMdatas)
    yesme <- function(gg) {length(which(gg == "Yes"))}
    
    ## checking where a subject belongs to
    
    wght <- do.call(rbind, lapply(seq_len(nrow(KMdatas)), function(ii) {
      
      data.frame("bad" = yesme(KMdatas[ii, colnames(KMdatas) %in% c(sa, ma)]), 
                 "neutral_bad" = yesme(KMdatas[ii, colnames(KMdatas) %in% c(na)]), 
                 "good" = yesme(KMdatas[ii, colnames(KMdatas) %in% c(sf, mf)]), 
                 "neutral_good" = yesme(KMdatas[ii, colnames(KMdatas) %in% c(nf)]))
    }))
    
    
    fav_ix <- which(wght$bad == 0 & wght$good > 0 )
    adv_ix <- which(wght$bad > 0 & wght$good == 0)
    int_ix <- which(wght$bad > 0 & wght$good > 0)
    int_fav_ix <- which(wght$good > 0 & wght$bad == 0 & wght$neutral_bad > 0)
    fav_ix <- which(wght$bad == 0 & wght$good > 0 )
    
    
    
    KMdata$MOD_RISK_PRACTICAL_V2 <- rep("Intermediate", nrow(KMdata))
    KMdata$MOD_RISK_PRACTICAL_V2[adv_ix] <- "Adverse"
    KMdata$MOD_RISK_PRACTICAL_V2[fav_ix] <- "Favorable"
    if(length(int_fav_ix) > 0){
      KMdata$MOD_RISK_PRACTICAL_V2[int_fav_ix] <- "Intermediate" # Favorable get trumped by Interemdiate due to FLT3ITD
    } 
    
    
    try(KMdata <- KMdata_test <- KMdata %>% mutate(MOD_RISK_PRACTICAL_V2 = fct_relevel(MOD_RISK_PRACTICAL_V2, c("Favorable", "Intermediate", "Adverse")), 
                                                   MOD_RISK_PRACTICAL_V2 = relevel(MOD_RISK_PRACTICAL_V2, ref = 'Favorable')))
    KMdata_sct$MOD_RISK_PRACTICAL_V2 <- KMdata$MOD_RISK_PRACTICAL_V2
    
    
    
    ## 2nd def
    sa_edit <- sa_edit[sa_edit %in% c(grep("CYT_", sa_edit, value = T), grep("TP53_", sa_edit, value = T))]
    ma_edit <- ma_edit[ma_edit %in% c(grep("CYT_", ma_edit, value = T), grep("TP53_", ma_edit, value = T))] 
    
    sf_edit <- sf_edit[sf_edit %in% c(grep("CYT_", sf_edit, value = T), grep("IDH1", sf_edit, value = T), grep("IDH2", sf_edit, value = T), grep("NPM1", sf_edit, value = T))] 
    mf_edit <- mf_edit[mf_edit %in% c(grep("CYT_", mf_edit, value = T), grep("IDH1", mf_edit, value = T), grep("IDH2", mf_edit, value = T), grep("NPM1", mf_edit, value = T))] 
    
    na_edit <- "FLT3ITD_MUTATION" #c("NGS_FLT3ITD", "PCR_FLT3_ITD")
    
    KMdatas <- KMdata[, colnames(KMdata) %in% c(sa_edit, ma_edit, mf_edit, sf_edit, na_edit, nf_edit, "time", "censor", "ELN_RISK_GROUP",
                                                "TP53_MUTATION", "NGS_KRAS", "NGS_NRAS")]
    
    dim(KMdatas)
    yesme <- function(gg) {length(which(gg == "Yes"))}
    
    ## checking where a subject belongs to
    
    wght_edit <- do.call(rbind, lapply(seq_len(nrow(KMdatas)), function(ii) {
      data.frame("bad" = yesme(KMdatas[ii, colnames(KMdatas) %in% c(sa_edit, ma_edit)]), 
                 "neutral_bad" = yesme(KMdatas[ii, colnames(KMdatas) %in% c(na_edit)]), 
                 "good" = yesme(KMdatas[ii, colnames(KMdatas) %in% c(sf_edit, mf_edit)]), 
                 "neutral_good" = yesme(KMdatas[ii, colnames(KMdatas) %in% c(nf_edit)]))
    }))
    
    
    fav_edit_ix <- which(wght_edit$bad == 0 & wght_edit$good > 0 )
    adv_edit_ix <- which(wght_edit$bad > 0 & wght_edit$good == 0)
    int_edit_ix <- which(wght_edit$bad > 0 & wght_edit$good > 0)
    int_fav_edit_ix <- which(wght_edit$good > 0 & wght_edit$bad == 0 & wght_edit$neutral_bad > 0)
    fav_edit_ix <- which(wght_edit$bad == 0 & wght_edit$good > 0 )
    
    
    
    KMdata$MOD_RISK_EDIT_PRACTICAL_V2 <- rep("Intermediate", nrow(KMdata))
    KMdata$MOD_RISK_EDIT_PRACTICAL_V2[adv_edit_ix] <- "Adverse"
    KMdata$MOD_RISK_EDIT_PRACTICAL_V2[fav_edit_ix] <- "Favorable"
    if(length(int_fav_edit_ix) > 0){
      KMdata$MOD_RISK_EDIT_PRACTICAL_V2[int_fav_edit_ix] <- "Intermediate" # Favorable get trumped by Interemdiate due to FLT3ITD
    } 
    
    
    try(KMdata <- KMdata_test <- KMdata %>% mutate(MOD_RISK_EDIT_PRACTICAL_V2 = fct_relevel(MOD_RISK_EDIT_PRACTICAL_V2, c("Favorable", "Intermediate", "Adverse")), 
                                                   MOD_RISK_EDIT_PRACTICAL_V2 = relevel(MOD_RISK_EDIT_PRACTICAL_V2, ref = 'Favorable')))
    KMdata_sct$MOD_RISK_EDIT_PRACTICAL_V2 <- KMdata$MOD_RISK_EDIT_PRACTICAL_V2
  }
  
  

  col_set <- col_set2 <- c("green", "maroon", "blue")
  
  w1 <- 650; h1 = 430
  
  
  margin <- 0.30; 
  minn <- 5; 
  
  nw_nm <- rep("Risk group", 2)
  
  if(is.null(global_var) == FALSE) {
    nw <- c("MOD_RISK_PRACTICAL_V2", "MOD_RISK_EDIT_PRACTICAL_V2") # used to name plots before saving files
    KMdata$MOD_RISK <- KMdata$MOD_RISK_PRACTICAL_V2
    KMdata_sct$MOD_RISK <- KMdata_sct$MOD_RISK_PRACTICAL_V2
    
    KMdata$MOD_RISK_EDIT <- KMdata$MOD_RISK_EDIT_PRACTICAL_V2
    KMdata_sct$MOD_RISK_EDIT <- KMdata_sct$MOD_RISK_EDIT_PRACTICAL_V2
    
  } else{
    nw <- c("MOD_RISK", "MOD_RISK_EDIT") 
  }
  
  
  table(KMdata$MOD_RISK)
  
 #data_type == "predictive (cv)" | data_type == "predictive (test)"
  
  if(is.null(test_ix) == FALSE & data_type == "predictive (cv)"){
    KMdata_test <- KMdata %>% filter(MRN %in% as.character(test_ix))
    KMdata_sct_test <- KMdata_sct %>% filter(MRN %in% as.character(test_ix))
    KMdata <- KMdata %>% filter(MRN %!in% as.character(test_ix))
    KMdata_sct <- KMdata_sct %>% filter(MRN %!in% as.character(test_ix))
    KMdata_sct$MOD_RISK <- KMdata$MOD_RISK
    KMdata_sct$MOD_RISK_EDIT <- KMdata$MOD_RISK_EDIT
  } 
  
   
    risk_df <- test_risk <- test_mrn <- res_surv <- t_risksetROC <- t_survivalROC <- ibrier <-
      t_risksetROC_eln <- t_survivalROC_eln <- ibrier_eln <- brier_365_eln <-  brier_365 <-
     
      ibrier_cridge <- t_survivalROC_ridge <- t_risksetROC_ridge <-
      ibrier_cridge_eln <-  t_survivalROC_ridge_eln <-   t_risksetROC_ridge_eln <- NULL

      risk_df_edit <- test_risk_edit <- res_surv_eln <- res_surv <-
      res_surv_edit <- t_risksetROC_edit <- t_survivalROC_edit <- ibrier_edit <-  
      brier_365_edit <- ibrier_cridge_edit <- t_survivalROC_ridge_edit <- t_risksetROC_ridge_edit <- NULL
  
  try(rm(jj))
  
  for(jj in 1 : length(nw)){ 
    
    icl <- NA; icl_mrn <- NULL
    
    
    # complete cases 
    if( length(grep("complete", os_type, value = T)) > 0 ){
      
      if(jj == 1){
        ll <- c(sf, mf, sa, ma) 
      } else if (jj == 2) {
        ll <- c(sf_edit, mf_edit, sa_edit, ma_edit) 
      }
      ll <- ll[which(ll %in% colnames(KMdata) == TRUE)]
      icl <- paste0("which(" , paste0("KMdata$", ll, " == 'Missing' ", collapse = " | "), ")")
      icl <- eval(parse(text = icl))
      icl_mrn <- as.character(KMdata$MRN[icl])
      
    }
    
    
    
    ## ITT : all pts
    try(run <- try(plotInput2_paperML( 
                              conf_int = TRUE, median_line = TRUE, 
                              test1 = "Equality", 
                              test2 = c("Log-rank", "Tarone-Ware", "Fleming-Harrington", "Renyi", "Max-Combo",  "RMST", "KONP", "mdir2", "mdir4"),
                              test_p = c(0, 0, 1), 
                              test_q = c(0, 0, 0),
                              col_set2, 
                              type_surv_prob = "Overall survival (OS)", 
                              name_var = nw_nm[jj], 
                            
                              xsurv = nw[jj],   
                              xsurv2 = "none", 
                              xsurv3 = "none",
                              Main = KMdata, 
                             
                              font.ytickslab = 17, 
                              font.xtickslab = 17,
                              font.title = 1,
                              font.x = 15, 
                              font.y = 15, 
                              censor.size = 5,
                              xlab = "Time (days) post-treatment",
                              conf.int.alpha = 0.2, 
                              conf.int.style = "ribbon",
                              legend_name = NULL, 
                              wghts = NULL)))  
    
    
    if(length(grep("complete", os_type, value = T)) > 0) {
    try(run_cc <- try(plotInput2_paperML( 
                              conf_int = TRUE, median_line = TRUE, 
                              test1 = "Equality", 
                              test2 = c("Log-rank", "Tarone-Ware", "Fleming-Harrington", "Renyi", "Max-Combo",  "RMST", "KONP", "mdir2", "mdir4"),
                              test_p = c(0, 0, 1), 
                              test_q = c(0, 0, 0),
                              col_set2, 
                              type_surv_prob = "Overall survival (OS)", 
                              name_var = nw_nm[jj], 
                              
                              xsurv = nw[jj],   
                              xsurv2 = "none", 
                              xsurv3 = "none",
                              Main = KMdata[seq_len(nrow(KMdata)) %!in% icl, ], 
                             
                              font.ytickslab = 17, 
                              font.xtickslab = 17,
                              font.title = 1,
                              font.x = 15, 
                              font.y = 15, 
                              censor.size = 5,
                              xlab = "Time (days) post-treatment",
                              conf.int.alpha = 0.2, 
                              conf.int.style = "ribbon",
                              legend_name = NULL, 
                              wght = NULL)))  
    }
    
    
    
    ## excluding sct 
    try(run_exsct <- try(plotInput2_paperML(censor_time = TRUE, over_day0 = "Treatment (first induction) start date", 
                                    over_day_cens = as.Date("2025-01-22"), 
                                    conf_int = TRUE, median_line = TRUE, 
                                    test1 = "Equality", 
                                    test2 = c("Log-rank", "Tarone-Ware", "Fleming-Harrington", "Renyi", "Max-Combo",  "RMST", "KONP", "mdir2", "mdir4"),
                                    test_p = c(0, 0, 1), 
                                    test_q = c(0, 0, 0),
                                    col_set2, 
                                    type_surv_prob = "Overall survival (OS)", 
                                    name_var = nw_nm[jj], 
                                   
                                    xsurv = nw[jj],   
                                    xsurv2 = "none", 
                                    xsurv3 = "none",
                                    Main = KMdata %>% subset(SCT_FLAG %in% c("ven/aza")), 
                                    
                                    font.ytickslab = 17, 
                                    font.xtickslab = 17,
                                    font.title = 1,
                                    font.x = 15, 
                                    font.y = 15, 
                                    censor.size = 5,
                                    xlab = "Time (days) post-treatment",
                                    conf.int.alpha = 0.2, 
                                    conf.int.style = "ribbon",
                                    legend_name = NULL, 
                                    wght = NULL)))  
    
    
    if(length(grep("complete", os_type, value = T)) > 0){
      try(run_cc_exsct <- try(plotInput2_paperML( 
                                      conf_int = TRUE, median_line = TRUE, 
                                      test1 = "Equality", 
                                      test2 = c("Log-rank", "Tarone-Ware", "Fleming-Harrington", "Renyi", "Max-Combo",  "RMST", "KONP", "mdir2", "mdir4"),
                                      test_p = c(0, 0, 1), 
                                      test_q = c(0, 0, 0),
                                      col_set2,  
                                      type_surv_prob = "Overall survival (OS)", 
                                      name_var = nw_nm[jj], 
                                      
                                      xsurv = nw[jj],   
                                      xsurv2 = "none", 
                                      xsurv3 = "none",
                                      Main = KMdata[seq_len(nrow(KMdata)) %!in% icl, ] %>% subset(SCT_FLAG %in% c("ven/aza")), 
                                     
                                      font.ytickslab = 17, 
                                      font.xtickslab = 17,
                                      font.title = 1,
                                      font.x = 15, 
                                      font.y = 15, 
                                      censor.size = 5,
                                      xlab = "Time (days) post-treatment",
                                      conf.int.alpha = 0.2, 
                                      conf.int.style = "ribbon",
                                      legend_name = NULL, 
                                      wght = NULL))) 
    }
    
    ## just including sct
    try(run_insct <- try(plotInput2_paperML( 
                                    conf_int = TRUE, median_line = TRUE, 
                                    test1 = "Equality", 
                                    test2 = c("Log-rank", "Tarone-Ware", "Fleming-Harrington", "Renyi", "Max-Combo",  "RMST", "KONP", "mdir2", "mdir4"),
                                    test_p = c(0, 0, 1), 
                                    test_q = c(0, 0, 0),
                                    col_set2, 
                                    type_surv_prob = "Overall survival (OS)", 
                                    name_var = nw_nm[jj], 
                                    
                                    xsurv = nw[jj],   
                                    xsurv2 = "none", 
                                    xsurv3 = "none",
                                    Main = KMdata %>% subset(SCT_FLAG %in% c("ven/aza + SCT")), 
                                  
                                    font.ytickslab = 17, 
                                    font.xtickslab = 17,
                                    font.title = 1,
                                    font.x = 15, 
                                    font.y = 15, 
                                    censor.size = 5,
                                    xlab = "Time (days) post-treatment",
                                    conf.int.alpha = 0.2, 
                                    conf.int.style = "ribbon",
                                    legend_name = NULL, 
                                    wght = NULL))) 
    
    if(length(grep("complete", os_type, value = T)) > 0){
      try(run_cc_insct <- try(plotInput2_paperML( 
                                      conf_int = TRUE, median_line = TRUE, 
                                      test1 = "Equality", 
                                      test2 = c("Log-rank", "Tarone-Ware", "Fleming-Harrington", "Renyi", "Max-Combo",  "RMST", "KONP", "mdir2", "mdir4"),
                                      test_p = c(0, 0, 1), 
                                      test_q = c(0, 0, 0),
                                      col_set2, 
                                      type_surv_prob = "Overall survival (OS)", 
                                      name_var = nw_nm[jj], 
                                     
                                      xsurv = nw[jj],   
                                      xsurv2 = "none", 
                                      xsurv3 = "none",
                                      Main = KMdata %>% subset(SCT_FLAG %in% c("ven/aza + SCT")), 
                                     
                                      font.ytickslab = 17, 
                                      font.xtickslab = 17,
                                      font.title = 1,
                                      font.x = 15, 
                                      font.y = 15, 
                                      censor.size = 5,
                                      xlab = "Time (days) post-treatment",
                                      conf.int.alpha = 0.2, 
                                      conf.int.style = "ribbon",
                                      legend_name = NULL, 
                                      wght = NULL))) 
      
      
    }
    
    ## treating sct as censored 
    try(run_sct <- try(plotInput2_paperML( 
                                  conf_int = TRUE, median_line = TRUE, 
                                  test1 = "Equality", 
                                  test2 = c("Log-rank", "Tarone-Ware", "Fleming-Harrington", "Renyi", "Max-Combo",  "RMST", "KONP", "mdir2", "mdir4"),
                                  test_p = c(0, 0, 1), 
                                  test_q = c(0, 0, 0),
                                  col_set2,      
                                  type_surv_prob = "Overall survival (OS)", 
                                  name_var = nw_nm[jj], 
                                 
                                  xsurv = nw[jj],   
                                  xsurv2 = "none", 
                                  xsurv3 = "none",
                                  Main = KMdata_sct %>% mutate(time == time_sct, censor == censor_sct), 
                               
                                  font.ytickslab = 17, 
                                  font.xtickslab = 17,
                                  font.title = 1,
                                  font.x = 15, 
                                  font.y = 15, 
                                  censor.size = 5,
                                  xlab = "Time (days) post-treatment",
                                  conf.int.alpha = 0.2, 
                                  conf.int.style = "ribbon",
                                  legend_name = NULL, 
                                  wght = NULL)))  
    

    if(length(grep("complete", os_type, value = T)) > 0){
      try(run_cc_sct <- try(plotInput2_paperML(censor_time = TRUE, over_day0 = "Treatment (first induction) start date", 
                                            over_day_cens = as.Date("2025-01-22"), 
                                            conf_int = TRUE, median_line = TRUE, 
                                            test1 = "Equality", 
                                            test2 = c("Log-rank", "Tarone-Ware", "Fleming-Harrington", "Renyi", "Max-Combo",  "RMST", "KONP", "mdir2", "mdir4"),
                                            test_p = c(0, 0, 1), 
                                            test_q = c(0, 0, 0),
                                            col_set2,      
                                            type_surv_prob = "Overall survival (OS)", 
                                            name_var = nw_nm[jj], 
                                            
                                            xsurv = nw[jj],   
                                            xsurv2 = "none", 
                                            xsurv3 = "none",
                                            Main = KMdata_sct[seq_len(nrow(KMdata_sct)) %!in% icl, ] %>% mutate(time == time_sct, censor == censor_sct), 
                                            
                                            font.ytickslab = 17, 
                                            font.xtickslab = 17,
                                            font.title = 1,
                                            font.x = 15, 
                                            font.y = 15, 
                                            censor.size = 5,
                                            xlab = "Time (days) post-treatment",
                                            conf.int.alpha = 0.2, 
                                            conf.int.style = "ribbon",
                                            legend_name = NULL, 
                                            wght = NULL)))  
    }
    
   
    # a1 <- run$prnt$g1$plot; b1 <- run$prnt$g1$table; a2 <- run_exsct$prnt$g1$plot; b2 <- run_exsct$prnt$g1$table; 
    # a3 <- run_sct$prnt$g1$plot; b3 <- run_sct$prnt$g1$table;
    # g1 <- run$prnt$g3; 
    
    
    a1 <- a2 <- a3 <- a4 <- a5 <-  b1 <- b2 <- b3 <- b4 <- b5 <- c1 <- c2 <- c3 <- c4 <- c5 <- d1 <- d2 <- d3 <- d4 <- d5 <- NULL
    if(class(run) != "try-error" & class(run_exsct) != "try-error" & class(run_sct) != "try-error") {
      g1 <- run$prnt$g3; t1 <- run$prnt$g1$table$data; a1 <- run$prnt$g1$plot; b1 <- run$prnt$g1$table; c1 <- c1_mod <- run$prnt$g1$data.survplot; 
      g2 <- run_exsct$prnt$g3; t2 <- run_exsct$prnt$g1$table$data; a2 <- run_exsct$prnt$g1$plot; b2 <- run_exsct$prnt$g1$table; c2 <- c2_mod <- run_exsct$prnt$g1$data.survplot; 
      g3 <- run_sct$prnt$g3; t3 <- run_sct$prnt$g1$table$data; a3 <- run_sct$prnt$g1$plot; b3 <- run_sct$prnt$g1$table; c3 <- c3_mod <- run_sct$prnt$g1$data.survplot;  a4 <- run_insct$prnt$g1$plot; b4 <- run_insct$prnt$g1$table; 
    }
    
    if(class(run) != "try-error" & class(run_exsct) == "try-error" & class(run_sct) != "try-error") {
      g1 <- run$prnt$g3; t1 <- run$prnt$g1$table$data; a1 <- run$prnt$g1$plot; b1 <- run$prnt$g1$table; c1 <- c1_mod <- run$prnt$g1$data.survplot; 
      g2 <- t2 <- a2 <- b2 <- c2 <- c2_mod <- ggplot() + theme_void() 
      g3 <- run_sct$prnt$g3; t3 <- run_sct$prnt$g1$table$data; a3 <- run_sct$prnt$g1$plot; b3 <- run_sct$prnt$g1$table; c3 <- c3_mod <- run_sct$prnt$g1$data.survplot; 
    }
    
    
    # now add the title
    library(cowplot)
    title <- ggdraw() + 
      draw_label(
        "",
        fontface = 'bold',
        x = 0,
        hjust = 0
      )
    png(paste0("", nw[jj],"_MERGED_FC_", fc_include,"_CD71_", exclude_cd71, "_", mis_type, "_", cen_type,"_", good_def, ".png"), width = 2700, height = 600)
    plot_list1 <- cowplot::plot_grid(plotlist = list(a1, a2, a3, a4,  
                                                     b1, b2, b3, b4), nrow = 2, ncol = 4, 
                                     byrow = T, rel_heights = c(6, 2), 
                                     label_size = 6, 
                                     rel_widths = c(4, 4, 4, 4), labels = c("All data", 
                                                                            "Excluding SCT", 
                                                                            "SCT as censored",
                                                                            "Including only SCT"))
    plot_list2 <- cowplot::plot_grid(title, plot_list1, nrow = 2, rel_heights = c(0.4, 2))
    plot_list2 <- cowplot::plot_grid(plot_list1, g1, ncol = 2, nrow = 1, rel_widths = c(5, 1))
    print(plot_list2)
    dev.off()
    
    
    if(length(grep("complete", os_type, value = T)) > 0) {
      
      g1_cc <- run_cc$prnt$g3; t1_cc <- run_cc$prnt$g1$table$data; a1_cc <- run_cc$prnt$g1$plot; b1_cc <- run_cc$prnt$g1$table; c1_cc <- run_cc$prnt$g1$data.survplot; 
      g2_cc <- run_cc_exsct$prnt$g3; t2_cc <- run_cc_exsct$prnt$g1$table$data; a2_cc <- run_cc_exsct$prnt$g1$plot; b2_cc <- run_cc_exsct$prnt$g1$table; c2_cc <- run_cc_exsct$prnt$g1$data.survplot; 
      g3_cc <- run_cc_sct$prnt$g3; t3_cc <- run_cc_sct$prnt$g1$table$data; a3_cc <- run_cc_sct$prnt$g1$plot; b3_cc <- run_cc_sct$prnt$g1$table; c3_cc <- run_cc_sct$prnt$g1$data.survplot; 
      g4_cc <- run_cc_insct$prnt$g3; t4_cc <- run_cc_insct$prnt$g1$table$data; a4_cc <- run_cc_insct$prnt$g1$plot; b4_cc <- run_cc_insct$prnt$g1$table; c4_cc <- run_cc_insct$prnt$g1$data.survplot; 
      
      png(paste0("COMPLETE_", nw[jj],"_MERGED_FC_", fc_include,"_CD71_", exclude_cd71, "_", mis_type, "_", cen_type,"_", good_def, ".png"), width = 2700, height = 600)
      plot_list1 <- cowplot::plot_grid(plotlist = list(a1_cc, a2_cc, a3_cc, a4_cc,  
                                                       b1_cc, b2_cc, b3_cc, b4_cc), nrow = 2, ncol = 4, 
                                       byrow = T, rel_heights = c(6, 2), 
                                       label_size = 5, 
                                       rel_widths = c(4, 4, 4, 4), labels = c("All data", 
                                                                              "Excluding SCT", 
                                                                              "SCT as censored",
                                                                              "Including only SCT"))
      plot_list2 <- cowplot::plot_grid(title, plot_list1, nrow = 2, rel_heights = c(0.4, 2))
      plot_list2 <- cowplot::plot_grid(plot_list1, g1, ncol = 2, nrow = 1, rel_widths = c(5, 1))
      print(plot_list2)
      dev.off()
    }
    
    
    ## ELN 
    if(length(grep("ELN", os_type, value = T)) > 0) {
    try(run_eln <- try(plotInput2_paperML( 
                                          conf_int = TRUE, median_line = TRUE, 
                                          test1 = "Equality", 
                                          test2 = c("Log-rank", "Tarone-Ware", "Fleming-Harrington", "Renyi", "Max-Combo",  "RMST", "KONP", "mdir2", "mdir4"),
                                          test_p = c(0, 0, 1), 
                                          test_q = c(0, 0, 0),
                                          col_set2, 
                                          type_surv_prob = "Overall survival (OS)", 
                                          name_var = "ELN", 
                                          
                                          xsurv = "ELN_RISK_GROUP",   
                                          xsurv2 = "none", 
                                          xsurv3 = "none",
                                          Main = KMdata, 
                                          
                                          font.ytickslab = 17, 
                                          font.xtickslab = 17,
                                          font.title = 1,
                                          font.x = 15, 
                                          font.y = 15, 
                                          censor.size = 5,
                                          xlab = "Time (days) post-treatment",
                                          conf.int.alpha = 0.2, 
                                          conf.int.style = "ribbon",
                                          legend_name = NULL, 
                                          wghts = NULL)))  
    
    if(length(grep("complete", os_type, value = T)) > 0) {
      try(run_eln_cc <- try(plotInput2_paperML( 
                                               conf_int = TRUE, median_line = TRUE, 
                                               test1 = "Equality", 
                                               test2 = c("Log-rank", "Tarone-Ware", "Fleming-Harrington", "Renyi", "Max-Combo",  "RMST", "KONP", "mdir2", "mdir4"),
                                               test_p = c(0, 0, 1), 
                                               test_q = c(0, 0, 0),
                                               col_set2, 
                                               type_surv_prob = "Overall survival (OS)", 
                                               name_var = "ELN", 
                                               
                                               xsurv = "ELN_RISK_GROUP",   
                                               xsurv2 = "none", 
                                               xsurv3 = "none",
                                               Main = KMdata %>% filter(ELN_RISK_GROUP != "Missing"), 
                                               
                                               font.ytickslab = 17, 
                                               font.xtickslab = 17,
                                               font.title = 1,
                                               font.x = 15, 
                                               font.y = 15, 
                                               censor.size = 5,
                                               xlab = "Time (days) post-treatment",
                                               conf.int.alpha = 0.2, 
                                               conf.int.style = "ribbon",
                                               legend_name = NULL, 
                                               wght = NULL)))  
    }
    
    
    
    ## excluding sct 
    try(run_eln_exsct <- try(plotInput2_paperML( 
                                                conf_int = TRUE, median_line = TRUE, 
                                                test1 = "Equality", 
                                                test2 = c("Log-rank", "Tarone-Ware", "Fleming-Harrington", "Renyi", "Max-Combo",  "RMST", "KONP", "mdir2", "mdir4"),
                                                test_p = c(0, 0, 1), 
                                                test_q = c(0, 0, 0),
                                                col_set2, 
                                                type_surv_prob = "Overall survival (OS)", 
                                                name_var = "ELN", 
                                                
                                                xsurv = "ELN_RISK_GROUP",   
                                                xsurv2 = "none", 
                                                xsurv3 = "none",
                                                Main = KMdata %>% subset(SCT_FLAG %in% c("ven/aza")), 
                                                
                                                font.ytickslab = 17, 
                                                font.xtickslab = 17,
                                                font.title = 1,
                                                font.x = 15, 
                                                font.y = 15, 
                                                censor.size = 5,
                                                xlab = "Time (days) post-treatment",
                                                conf.int.alpha = 0.2, 
                                                conf.int.style = "ribbon",
                                                legend_name = NULL, 
                                                wght = NULL)))  
    
    
    if(length(grep("complete", os_type, value = T)) > 0){
      try(run_eln_cc_exsct <- try(plotInput2_paperML( 
                                                     conf_int = TRUE, median_line = TRUE, 
                                                     test1 = "Equality", 
                                                     test2 = c("Log-rank", "Tarone-Ware", "Fleming-Harrington", "Renyi", "Max-Combo",  "RMST", "KONP", "mdir2", "mdir4"),
                                                     test_p = c(0, 0, 1), 
                                                     test_q = c(0, 0, 0),
                                                     col_set2,  
                                                     type_surv_prob = "Overall survival (OS)", 
                                                     name_var = "ELN", 
                                                     
                                                     xsurv = "ELN_RISK_GROUP",   
                                                     xsurv2 = "none", 
                                                     xsurv3 = "none",
                                                     Main = KMdata %>% subset(SCT_FLAG %in% c("ven/aza")) %>% filter(ELN_RISK_GROUP != "Missing"), 
                                                     
                                                     font.ytickslab = 17, 
                                                     font.xtickslab = 17,
                                                     font.title = 1,
                                                     font.x = 15, 
                                                     font.y = 15, 
                                                     censor.size = 5,
                                                     xlab = "Time (days) post-treatment",
                                                     conf.int.alpha = 0.2, 
                                                     conf.int.style = "ribbon",
                                                     legend_name = NULL, 
                                                     wght = NULL))) 
    }
    
    ## just including sct
    try(run_eln_insct <- try(plotInput2_paperML( 
                                                conf_int = TRUE, median_line = TRUE, 
                                                test1 = "Equality", 
                                                test2 = c("Log-rank", "Tarone-Ware", "Fleming-Harrington", "Renyi", "Max-Combo",  "RMST", "KONP", "mdir2", "mdir4"),
                                                test_p = c(0, 0, 1), 
                                                test_q = c(0, 0, 0),
                                                col_set2, 
                                                type_surv_prob = "Overall survival (OS)", 
                                                name_var = "ELN", 
                                                
                                                xsurv = "ELN_RISK_GROUP",   
                                                xsurv2 = "none", 
                                                xsurv3 = "none",
                                                Main = KMdata %>% subset(SCT_FLAG %in% c("ven/aza + SCT")), 
                                                
                                                font.ytickslab = 17, 
                                                font.xtickslab = 17,
                                                font.title = 1,
                                                font.x = 15, 
                                                font.y = 15, 
                                                censor.size = 5,
                                                xlab = "Time (days) post-treatment",
                                                conf.int.alpha = 0.2, 
                                                conf.int.style = "ribbon",
                                                legend_name = NULL, 
                                                wght = NULL))) 
    
    if(length(grep("complete", os_type, value = T)) > 0){
      try(run_eln_cc_insct <- try(plotInput2_paperML( 
                                                     conf_int = TRUE, median_line = TRUE, 
                                                     test1 = "Equality", 
                                                     test2 = c("Log-rank", "Tarone-Ware", "Fleming-Harrington", "Renyi", "Max-Combo",  "RMST", "KONP", "mdir2", "mdir4"),
                                                     test_p = c(0, 0, 1), 
                                                     test_q = c(0, 0, 0),
                                                     col_set2, 
                                                     type_surv_prob = "Overall survival (OS)", 
                                                     name_var = "ELN", 
                                                     
                                                     xsurv = "ELN_RISK_GROUP",   
                                                     xsurv2 = "none", 
                                                     xsurv3 = "none",
                                                     Main = KMdata %>% subset(SCT_FLAG %in% c("ven/aza + SCT")) %>% filter(ELN_RISK_GROUP != "Missing"), 
                                                     
                                                     font.ytickslab = 17, 
                                                     font.xtickslab = 17,
                                                     font.title = 1,
                                                     font.x = 15, 
                                                     font.y = 15, 
                                                     censor.size = 5,
                                                     xlab = "Time (days) post-treatment",
                                                     conf.int.alpha = 0.2, 
                                                     conf.int.style = "ribbon",
                                                     legend_name = NULL, 
                                                     wght = NULL))) 
      
      
    }
    
    ## treating sct as censored 
    try(run_eln_sct <- try(plotInput2_paperML( 
                                              conf_int = TRUE, median_line = TRUE, 
                                              test1 = "Equality", 
                                              test2 = c("Log-rank", "Tarone-Ware", "Fleming-Harrington", "Renyi", "Max-Combo",  "RMST", "KONP", "mdir2", "mdir4"),
                                              test_p = c(0, 0, 1), 
                                              test_q = c(0, 0, 0),
                                              col_set2,      
                                              type_surv_prob = "Overall survival (OS)", 
                                              name_var = "ELN", 
                                              
                                              xsurv = "ELN_RISK_GROUP",   
                                              xsurv2 = "none", 
                                              xsurv3 = "none",
                                              Main = KMdata_sct %>% mutate(time == time_sct, censor == censor_sct), 
                                              
                                              font.ytickslab = 17, 
                                              font.xtickslab = 17,
                                              font.title = 1,
                                              font.x = 15, 
                                              font.y = 15, 
                                              censor.size = 5,
                                              xlab = "Time (days) post-treatment",
                                              conf.int.alpha = 0.2, 
                                              conf.int.style = "ribbon",
                                              legend_name = NULL, 
                                              wght = NULL)))  
    
    
    if(length(grep("complete", os_type, value = T)) > 0){
      try(run_eln_cc_sct <- try(plotInput2_paperML( 
                                                   conf_int = TRUE, median_line = TRUE, 
                                                   test1 = "Equality", 
                                                   test2 = c("Log-rank", "Tarone-Ware", "Fleming-Harrington", "Renyi", "Max-Combo",  "RMST", "KONP", "mdir2", "mdir4"),
                                                   test_p = c(0, 0, 1), 
                                                   test_q = c(0, 0, 0),
                                                   col_set2,      
                                                   type_surv_prob = "Overall survival (OS)", 
                                                   name_var = "ELN", 
                                                   
                                                   xsurv = "ELN_RISK_GROUP",   
                                                   xsurv2 = "none", 
                                                   xsurv3 = "none",
                                                   Main = KMdata_sct[seq_len(nrow(KMdata_sct)) %!in% icl, ] %>% mutate(time == time_sct, censor == censor_sct) %>% filter(ELN_RISK_GROUP != "Missing"), 
                                                   
                                                   font.ytickslab = 17, 
                                                   font.xtickslab = 17,
                                                   font.title = 1,
                                                   font.x = 15, 
                                                   font.y = 15, 
                                                   censor.size = 5,
                                                   xlab = "Time (days) post-treatment",
                                                   conf.int.alpha = 0.2, 
                                                   conf.int.style = "ribbon",
                                                   legend_name = NULL, 
                                                   wght = NULL)))  
    }
    
    
    ## Printing ELN results
    png(paste0("ELN_MERGED_FC_", fc_include,"_CD71_", exclude_cd71, "_", mis_type, "_", cen_type,"_", good_def, ".png"), width = 2700, height = 600)
    g1 <- run_eln$prnt$g3; 
    c1_eln <- run_eln$prnt$g1$data.survplot
    c2_eln <- run_eln_exsct$prnt$g1$data.survplot
    c3_eln <- run_eln_sct$prnt$g1$data.survplot
    
    plot_list1 <- cowplot::plot_grid(plotlist = list(run_eln$prnt$g1$plot, run_eln_exsct$prnt$g1$plot, run_eln_sct$prnt$g1$plot, run_eln_insct$prnt$g1$plot,  
                                                     run_eln$prnt$g1$table, run_eln_exsct$prnt$g1$table,  run_eln_sct$prnt$g1$table, run_eln_insct$prnt$g1$table), nrow = 2, ncol = 4, 
                                     byrow = T, rel_heights = c(6, 2), 
                                     label_size = 5, 
                                     rel_widths = c(4, 4, 4, 4), labels = c("All data", 
                                                                            "Excluding SCT", 
                                                                            "SCT as censored",
                                                                            "Including only SCT"))
    plot_list2 <- cowplot::plot_grid(title, plot_list1, nrow = 2, rel_heights = c(0.4, 2))
    plot_list2 <- cowplot::plot_grid(plot_list1, g1, ncol = 2, nrow = 1, rel_widths = c(5, 1))
    print(plot_list2)
    dev.off()
    
    
    ## ELN Complete 
    if(length(grep("complete", os_type, value = T)) > 0){
    png(paste0("ELN_COMPLETE_MERGED_FC_", fc_include,"_CD71_", exclude_cd71, "_", mis_type, "_", cen_type,"_", good_def, ".png"), width = 2700, height = 600)
    plot_list1 <- cowplot::plot_grid(plotlist = list(run_eln_cc$prnt$g1$plot, run_eln_cc_exsct$prnt$g1$plot, run_eln_cc_sct$prnt$g1$plot, run_eln_cc_insct$prnt$g1$plot,  
                                                       run_eln_cc$prnt$g1$table, run_eln_cc_exsct$prnt$g1$table, run_eln_cc_sct$prnt$g1$table, run_eln_cc_insct$prnt$g1$table), nrow = 2, ncol = 4, 
                                       byrow = T, rel_heights = c(6, 2), 
                                       label_size = 5, 
                                       rel_widths = c(4, 4, 4, 4), labels = c("All data", 
                                                                              "Excluding SCT", 
                                                                              "SCT as censored",
                                                                              "Including only SCT"))
      plot_list2 <- cowplot::plot_grid(title, plot_list1, nrow = 2, rel_heights = c(0.4, 2))
      plot_list2 <- cowplot::plot_grid(plot_list1, g1, ncol = 2, nrow = 1, rel_widths = c(5, 1))
      print(plot_list2)
      dev.off()
    }
    } # if ELN needs to be estimated

    ## Pairwise comparisons iff both ELN and pairwise comparisons have been selected
    if(length(grep("pairwise", os_type, value = T)) > 0 & length(grep("ELN", os_type, value = T)) > 0){ # Do only when all 3 categories are identified
      
      KMdata_sct <- KMdata %>% mutate(time == time_sct, censor == censor_sct)
      xxn <- c("Adverse risk", "Intermediate risk", "Favorable risk")
      col_rsk <- list(c("khaki4", "salmon"),
                      c("darkred", "magenta"),
                      c("goldenrod4", "seagreen4"))
      
      # adverse
      if(jj == 1){         # method-1
        l1 <- which(KMdata$MOD_RISK == "Adverse")
      }  else if(jj == 2){ # method-2
        l1 <- which(KMdata$MOD_RISK_EDIT == "Adverse")   
      }
      
      l2 <- which(KMdata$ELN_RISK_GROUP == "Adverse")
      
      df1 <- data.frame("MRN" = c(KMdata$MRN[l1], KMdata$MRN[l2]),
                        "time" = c(KMdata$time[l1], KMdata$time[l2]),
                        "censor" = c(KMdata$censor[l1], KMdata$censor[l2]),
                        "Group"  = c(rep("Adverse-VA", length(l1)), rep("Adverse-ELN", length(l2))),
                        "SCT_FLAG" = c(as.character(KMdata$SCT_FLAG[l1]), as.character(KMdata$SCT_FLAG[l2])))
      
      # adverse _sct
      if(jj == 1 ){       # method-1
        l1 <- which(KMdata_sct$MOD_RISK == "Adverse")
      } else if(jj == 2){ # method-2
        l1 <- which(KMdata_sct$MOD_RISK_EDIT == "Adverse")
      }
      l2 <- which(KMdata_sct$ELN_RISK_GROUP == "Adverse")
      
      df1_sct <- data.frame("MRN" = c(KMdata_sct$MRN[l1], KMdata_sct$MRN[l2]),
                            "time" = c(KMdata_sct$time[l1], KMdata_sct$time[l2]),
                            "censor" = c(KMdata_sct$censor[l1], KMdata_sct$censor[l2]),
                            "Group"  = c(rep("Adverse-VA", length(l1)), rep("Adverse-ELN", length(l2))))
      
      
      # intermediate
      if(jj == 1){       # method-1
        l1 <- which(KMdata$MOD_RISK == "Intermediate")
      } else if(jj == 2){ # method-2
        l1 <- which(KMdata$MOD_RISK_EDIT == "Intermediate")
      }
      l2 <- which(KMdata$ELN_RISK_GROUP == "Intermediate")
      
      df2 <- data.frame("MRN" = c(KMdata$MRN[l1], KMdata$MRN[l2]),
                        "time" = c(KMdata$time[l1], KMdata$time[l2]),
                        "censor" = c(KMdata$censor[l1], KMdata$censor[l2]),
                        "Group"  = c(rep("Intermediate-VA", length(l1)), rep("Intermediate-ELN", length(l2))),
                        "SCT_FLAG" = c(as.character(KMdata$SCT_FLAG[l1]), as.character(KMdata$SCT_FLAG[l2])))
      
      
      # intermediate_sct
      if(jj == 1){        # method-1
        l1 <- which(KMdata_sct$MOD_RISK == "Intermediate")
      } else if(jj == 2){ # method-2
        l1 <- which(KMdata_sct$MOD_RISK_EDIT == "Intermediate")
      }
      
      l2 <- which(KMdata_sct$ELN_RISK_GROUP == "Intermediate")
      df2_sct <- data.frame("MRN" = c(KMdata_sct$MRN[l1], KMdata_sct$MRN[l2]),
                            "time" = c(KMdata_sct$time[l1], KMdata_sct$time[l2]),
                            "censor" = c(KMdata_sct$censor[l1], KMdata_sct$censor[l2]),
                            "Group"  = c(rep("Intermediate-VA", length(l1)), rep("Intermediate-ELN", length(l2))))
      
      
      # Favorable
      if(jj == 1){        # method-1
        l1 <- which(KMdata$MOD_RISK == "Favorable")
      } else if(jj == 2){ # method-2
        l1 <- which(KMdata$MOD_RISK_EDIT == "Favorable")
      }
      
      l2 <- which(KMdata$ELN_RISK_GROUP == "Favorable")
      df3 <- data.frame("MRN" = c(KMdata$MRN[l1], KMdata$MRN[l2]),
                        "time" = c(KMdata$time[l1], KMdata$time[l2]),
                        "censor" = c(KMdata$censor[l1], KMdata$censor[l2]),
                        "Group"  = c(rep("Favorable-VA", length(l1)), rep("Favorable-ELN", length(l2))),
                        "SCT_FLAG" = c(as.character(KMdata$SCT_FLAG[l1]), as.character(KMdata$SCT_FLAG[l2])))
      
      
      # Favorable_sct
      if(jj == 1 ){        # method-1
        l1 <- which(KMdata_sct$MOD_RISK == "Favorable")
      } else if(jj == 2 ){ # method-2
        l1 <- which(KMdata_sct$MOD_RISK_EDIT == "Favorable")
      }
      l2 <- which(KMdata_sct$ELN_RISK_GROUP == "Favorable")
      df3_sct <- data.frame("MRN" = c(KMdata_sct$MRN[l1], KMdata_sct$MRN[l2]),
                            "time" = c(KMdata_sct$time[l1], KMdata_sct$time[l2]),
                            "censor" = c(KMdata_sct$censor[l1], KMdata_sct$censor[l2]),
                            "Group"  = c(rep("Favorable-VA", length(l1)), rep("Favorable-ELN", length(l2))))
      
      xxn <- c("Adverse risk", "Intermediate risk", "Favorable risk")
      col_rsk <- list(c("khaki4", "salmon"),
                      c("darkred", "magenta"),
                      c("goldenrod4", "seagreen4")
      )
      
      
      ## method-1
      for(xx in 1 : 3){
        
        QQ <- paste0("df <- df", xx, "")
        QQ <- eval(parse(text = QQ))
        
        QQ <- paste0("df_sct <- df", xx, "_sct")
        QQ <- eval(parse(text = QQ))
        
        
        try(run <- try(plotInput2_paperML(
                                          conf_int = TRUE, median_line = TRUE,
                                          test1 = "Equality",
                                          test2 = c("Log-rank", "Tarone-Ware", "Fleming-Harrington", "Renyi", "Max-Combo",  "RMST", "KONP", "mdir2", "mdir4"),
                                          test_p = c(0, 0, 1),
                                          test_q = c(0, 0, 0),
                                          col_rsk[[xx]],
                                          type_surv_prob = "Overall survival (OS)",
                                          name_var = xxn[xx],
                                          
                                          xsurv = "Group",
                                          xsurv2 = "none",
                                          xsurv3 = "none",
                                          Main = df,
                                          
                                          font.ytickslab = 17,
                                          font.xtickslab = 17,
                                          font.title = 1,
                                          font.x = 15,
                                          font.y = 15,
                                          censor.size = 5,
                                          xlab = "Time (days) post-treatment",
                                          conf.int.alpha = 0.2,
                                          conf.int.style = "ribbon",
                                          legend_name = NULL,
                                          wght = NULL))) 
        
        if(length(grep("complete", os_type, value = T)) > 0) { # original
          try(run_cc <- try(plotInput2_paperML(
                                               conf_int = TRUE, median_line = TRUE,
                                               test1 = "Equality",
                                               test2 = c("Log-rank", "Tarone-Ware", "Fleming-Harrington", "Renyi", "Max-Combo",  "RMST", "KONP", "mdir2", "mdir4"),
                                               test_p = c(0, 0, 1),
                                               test_q = c(0, 0, 0),
                                               col_rsk[[xx]],
                                               type_surv_prob = "Overall survival (OS)",
                                               name_var = xxn[xx],
                                               
                                               xsurv = "Group",
                                               xsurv2 = "none",
                                               xsurv3 = "none",
                                               Main = df %>% subset(MRN %!in% icl_mrn),
                                               
                                               font.ytickslab = 17,
                                               font.xtickslab = 17,
                                               font.title = 1,
                                               font.x = 15,
                                               font.y = 15,
                                               censor.size = 5,
                                               xlab = "Time (days) post-treatment",
                                               conf.int.alpha = 0.2,
                                               conf.int.style = "ribbon",
                                               legend_name = NULL,
                                               wght = NULL)))  
        }
        
        
        try(run_exsct <- try(plotInput2_paperML(
                                                conf_int = TRUE, median_line = TRUE,
                                                test1 = "Equality",
                                                test2 = c("Log-rank", "Tarone-Ware", "Fleming-Harrington", "Renyi", "Max-Combo",  "RMST", "KONP", "mdir2", "mdir4"),
                                                test_p = c(0, 0, 1),
                                                test_q = c(0, 0, 0),
                                                col_rsk[[xx]],
                                                type_surv_prob = "Overall survival (OS)",
                                                name_var = xxn[xx],
                                                
                                                xsurv = "Group",
                                                xsurv2 = "none",
                                                xsurv3 = "none",
                                                Main = df %>% subset(SCT_FLAG %in% c("ven/aza")),
                                                
                                                font.ytickslab = 17,
                                                font.xtickslab = 17,
                                                font.title = 1,
                                                font.x = 15,
                                                font.y = 15,
                                                censor.size = 5,
                                                xlab = "Time (days) post-treatment",
                                                conf.int.alpha = 0.2,
                                                conf.int.style = "ribbon",
                                                legend_name = NULL,
                                                wght = NULL)))  #for classic
        
        if(length(grep("complete", os_type, value = T)) > 0) { # original
          try(run_cc_exsct <- try(plotInput2_paperML(
                                                     conf_int = TRUE, median_line = TRUE,
                                                     test1 = "Equality",
                                                     test2 = c("Log-rank", "Tarone-Ware", "Fleming-Harrington", "Renyi", "Max-Combo",  "RMST", "KONP", "mdir2", "mdir4"),
                                                     test_p = c(0, 0, 1),
                                                     test_q = c(0, 0, 0),
                                                     col_rsk[[xx]],
                                                     type_surv_prob = "Overall survival (OS)",
                                                     name_var = xxn[xx],
                                                     
                                                     xsurv = "Group",
                                                     xsurv2 = "none",
                                                     xsurv3 = "none",
                                                     Main = df %>% subset(MRN %!in% icl_mrn) %>% subset(SCT_FLAG %in% c("ven/aza")),
                                                     
                                                     font.ytickslab = 17,
                                                     font.xtickslab = 17,
                                                     font.title = 1,
                                                     font.x = 15,
                                                     font.y = 15,
                                                     censor.size = 5,
                                                     xlab = "Time (days) post-treatment",
                                                     conf.int.alpha = 0.2,
                                                     conf.int.style = "ribbon",
                                                     legend_name = NULL,
                                                     wght = NULL))) 
        }
        
        try(run_sct <- try(plotInput2_paperML(
                                              conf_int = TRUE, median_line = TRUE,
                                              test1 = "Equality",
                                              test2 = c("Log-rank", "Tarone-Ware", "Fleming-Harrington", "Renyi", "Max-Combo",  "RMST", "KONP", "mdir2", "mdir4"),
                                              test_p = c(0, 0, 1),
                                              test_q = c(0, 0, 0),
                                              col_rsk[[xx]], 
                                              type_surv_prob = "Overall survival (OS)",
                                              name_var = xxn[[xx]],
                                              
                                              xsurv = "Group",
                                              xsurv2 = "none",
                                              xsurv3 = "none",
                                              Main = df_sct,
                                              
                                              font.ytickslab = 17,
                                              font.xtickslab = 17,
                                              font.title = 1,
                                              font.x = 15,
                                              font.y = 15,
                                              censor.size = 5,
                                              xlab = "Time (days) post-treatment",
                                              conf.int.alpha = 0.2,
                                              conf.int.style = "ribbon",
                                              legend_name = NULL,
                                              wght = NULL)))  
        
        
        if(length(grep("complete", os_type, value = T)) > 0) { 
          try(run_cc_sct <- try(plotInput2_paperML(
                                                   conf_int = TRUE, median_line = TRUE,
                                                   test1 = "Equality",
                                                   test2 = c("Log-rank", "Tarone-Ware", "Fleming-Harrington", "Renyi", "Max-Combo",  "RMST", "KONP", "mdir2", "mdir4"),
                                                   test_p = c(0, 0, 1),
                                                   test_q = c(0, 0, 0),
                                                   col_rsk[[xx]], 
                                                   type_surv_prob = "Overall survival (OS)",
                                                   name_var = xxn[[xx]],
                                                   
                                                   xsurv = "Group",
                                                   xsurv2 = "none",
                                                   xsurv3 = "none",
                                                   Main = df_sct %>% subset(MRN %!in% icl_mrn),
                                                   
                                                   font.ytickslab = 17,
                                                   font.xtickslab = 17,
                                                   font.title = 1,
                                                   font.x = 15,
                                                   font.y = 15,
                                                   censor.size = 5,
                                                   xlab = "Time (days) post-treatment",
                                                   conf.int.alpha = 0.2,
                                                   conf.int.style = "ribbon",
                                                   legend_name = NULL,
                                                   wght = NULL)))  #for classic
        }
        
        if(class(run) != "try-error" & class(run_exsct) != "try-error" & class(run_sct) != "try-error") {
          g1 <- run$prnt$g3; t1 <- run$prnt$g1$table$data; a1 <- run$prnt$g1$plot; b1 <- run$prnt$g1$table; c1 <- run$prnt$g1$data.survplot
          g2 <- run_exsct$prnt$g3; t2 <- run_exsct$prnt$g1$table$data; a2 <- run_exsct$prnt$g1$plot; b2 <- run_exsct$prnt$g1$table; c2 <- run_exsct$prnt$g1$data.survplot
          g3 <- run_sct$prnt$g3; t3 <- run_sct$prnt$g1$table$data; a3 <- run_sct$prnt$g1$plot; b3 <- run_sct$prnt$g1$table; c3 <- run_sct$prnt$g1$data.survplot
          
          
          png(paste0("", nw[jj],"_MERGED_", xxn[xx], "_FC_", fc_include,"_CD71_", exclude_cd71, "_", mis_type, "_", cen_type,"_", good_def, ".png"), width = 2000, height = 500)
          
          plot_list1 <- cowplot::plot_grid(plotlist = list(a1, a2, a3,
                                                           b1, b2, b3), nrow = 2, ncol = 3,
                                           byrow = T, rel_heights = c(6, 2),
                                           label_size = 6,
                                           rel_widths = c(4, 4, 4), labels = c("All data",
                                                                               "Excluding SCT",
                                                                               "SCT as censored"))
          
          plot_list2 <- cowplot::plot_grid(plot_list1, g1, ncol = 2, nrow = 1, rel_widths = c(5, 1))
          print(plot_list2)
          dev.off()
          
          if(jj == 1 ){ # After exlcuding any missing ELN
            ikp1 <- irr::kappam.fleiss(cbind(as.character(KMdata$MOD_RISK[which(KMdata$ELN_RISK_GROUP %!in% "Missing")]), as.character(KMdata$ELN_RISK_GROUP[which(KMdata$ELN_RISK_GROUP %!in% "Missing")])), exact = T)
            ikp2 <- irr::kappam.fleiss(cbind(as.character(KMdata$MOD_RISK[which(KMdata$ELN_RISK_GROUP %!in% "Missing")]), as.character(KMdata$ELN_RISK_GROUP[which(KMdata$ELN_RISK_GROUP %!in% "Missing")])), detail = T)
          }
          
          
          
          if(jj == 2 ){ # After exlcuding any missing ELN
            ikp1 <- irr::kappam.fleiss(cbind(as.character(KMdata$MOD_RISK_EDIT[which(KMdata$ELN_RISK_GROUP %!in% "Missing")]), as.character(KMdata$ELN_RISK_GROUP[which(KMdata$ELN_RISK_GROUP %!in% "Missing")])), exact = T)
            ikp2 <- irr::kappam.fleiss(cbind(as.character(KMdata$MOD_RISK_EDIT[which(KMdata$ELN_RISK_GROUP %!in% "Missing")]), as.character(KMdata$ELN_RISK_GROUP[which(KMdata$ELN_RISK_GROUP %!in% "Missing")])), detail = T)
          }
          
          saveRDS(list(ikp1 = ikp1, ikp2 = ikp2), file = paste0("", nw[jj],"_MERGED_KAPPA_FC_", fc_include,"_CD71_", exclude_cd71, "_", mis_type, "_", cen_type,"_", good_def, ".rds"))
          
        } 
        
        
        if(length(grep("complete", os_type, value = T)) > 0) {
          g1 <- run_cc$prnt$g3; t1 <- run_cc$prnt$g1$table$data; a1 <- run_cc$prnt$g1$plot; b1 <- run_cc$prnt$g1$table; c1 <- run_cc$prnt$g1$data.survplot
          g2 <- run_cc_exsct$prnt$g3; t2 <- run_cc_exsct$prnt$g1$table$data; a2 <- run_cc_exsct$prnt$g1$plot; b2 <- run_cc_exsct$prnt$g1$table; c2 <- run_cc_exsct$prnt$g1$data.survplot
          g3 <- run_cc_sct$prnt$g3; t3 <- run_cc_sct$prnt$g1$table$data; a3 <- run_cc_sct$prnt$g1$plot; b3 <- run_cc_sct$prnt$g1$table; c3 <- run_cc_sct$prnt$g1$data.survplot
          
          png(paste0("COMPLETE_", nw[jj],"_MERGED_", xxn[xx], "_FC_", fc_include,"_CD71_", exclude_cd71, "_", mis_type, "_", cen_type,"_", good_def, ".png"), width = 2000, height = 500)
          plot_list1 <- cowplot::plot_grid(plotlist = list(a1, a2, a3,
                                                           b1, b2, b3), nrow = 2, ncol = 3,
                                           byrow = T, rel_heights = c(6, 2),
                                           label_size = 6,
                                           rel_widths = c(4, 4, 4), labels = c("All data",
                                                                               "Excluding SCT",
                                                                               "SCT as censored"))
          plot_list2 <- cowplot::plot_grid(plot_list1, g1, ncol = 2, nrow = 1, rel_widths = c(5, 1))
          print(plot_list2)
          dev.off()
          
          
          if(jj == 1 ){
            ikp1 <- irr::kappam.fleiss(cbind(as.character(KMdata$MOD_RISK[which(KMdata$MRN %!in% icl_mrn | KMdata$ELN_RISK_GROUP %!in% "Missing")]), as.character(KMdata$ELN_RISK_GROUP[which(KMdata$MRN %!in% icl_mrn | KMdata$ELN_RISK_GROUP %!in% "Missing")])), exact = T)
            ikp2 <- irr::kappam.fleiss(cbind(as.character(KMdata$MOD_RISK[which(KMdata$MRN %!in% icl_mrn | KMdata$ELN_RISK_GROUP %!in% "Missing")]), as.character(KMdata$ELN_RISK_GROUP[which(KMdata$MRN %!in% icl_mrn | KMdata$ELN_RISK_GROUP %!in% "Missing")])), detail = T)
          }
          
          if(jj == 2 ){
            ikp1 <- irr::kappam.fleiss(cbind(as.character(KMdata$MOD_RISK_EDIT[which(KMdata$MRN %!in% icl_mrn | KMdata$ELN_RISK_GROUP %!in% "Missing")]), as.character(KMdata$ELN_RISK_GROUP[which(KMdata$MRN %!in% icl_mrn | KMdata$ELN_RISK_GROUP %!in% "Missing")])), exact = T)
            ikp2 <- irr::kappam.fleiss(cbind(as.character(KMdata$MOD_RISK_EDIT[which(KMdata$MRN %!in% icl_mrn | KMdata$ELN_RISK_GROUP %!in% "Missing")]), as.character(KMdata$ELN_RISK_GROUP[which(KMdata$MRN %!in% icl_mrn | KMdata$ELN_RISK_GROUP %!in% "Missing")])), detail = T)
          }
          
          saveRDS(list(ikp1 = ikp1, ikp2 = ikp2), file = paste0("COMPLETE_", nw[jj],"_MERGED_KAPPA_FC_", fc_include,"_CD71_", exclude_cd71, "_", mis_type, "_", cen_type,"_", good_def, ".rds"))
        }
        
        try(rm(a1)); try(rm(a2)); try(rm(a3));  try(rm(a4));
        try(rm(b1)); try(rm(b2)); try(rm(b3)); try(rm(b4));
        try(rm(g1)); try(rm(g2)); try(rm(g3)); try(rm(g4));
        try(rm(t1)); try(rm(t2)); try(rm(t3)); try(rm(t4));
      }
      
    }
    
    
    ## did not do it for complete cases within cv test dataset to save computing time BUT did it for imputed dataset
    if(is.null(test_ix) == FALSE & (data_type == "predictive (cv)" | data_type == "predictive (test)")){
      
      ## KM-specific prediction
      if(cen_type == "all"){ # ITT
        
        res_surv <- list(
          "Adverse" = data.frame("time" = c1_mod$time[which(c1_mod$MOD_RISK == "Adverse")], 
                                 "censor" = c1_mod$n.censor[which(c1_mod$MOD_RISK == "Adverse")], 
                                 "surv"  = c1_mod$surv[which(c1_mod$MOD_RISK == "Adverse")], 
                                 "lower" = c1_mod$lower[which(c1_mod$MOD_RISK == "Adverse")], 
                                 "upper" = c1_mod$upper[which(c1_mod$MOD_RISK == "Adverse")]), 
          
          "Favorable" = data.frame("time" = c1_mod$time[which(c1_mod$MOD_RISK == "Favorable")], 
                                   "censor" = c1_mod$n.censor[which(c1_mod$MOD_RISK == "Favorable")], 
                                   "surv"  = c1_mod$surv[which(c1_mod$MOD_RISK == "Favorable")], 
                                   "lower" = c1_mod$lower[which(c1_mod$MOD_RISK == "Favorable")], 
                                   "upper" = c1_mod$upper[which(c1_mod$MOD_RISK == "Favorable")]),
          
          "Intermediate" = data.frame("time" = c1_mod$time[which(c1_mod$MOD_RISK == "Intermediate")], 
                                      "censor" = c1_mod$n.censor[which(c1_mod$MOD_RISK == "Intermediate")], 
                                      "surv"  = c1_mod$surv[which(c1_mod$MOD_RISK == "Intermediate")], 
                                      "lower" = c1_mod$lower[which(c1_mod$MOD_RISK == "Intermediate")], 
                                      "upper" = c1_mod$upper[which(c1_mod$MOD_RISK == "Intermediate")])
        )
        
        res_surv_edit <- list(
          "Adverse" = data.frame("time" = c1_mod$time[which(c1_mod$MOD_RISK_EDIT == "Adverse")], 
                                 "censor" = c1_mod$n.censor[which(c1_mod$MOD_RISK_EDIT == "Adverse")], 
                                 "surv"  = c1_mod$surv[which(c1_mod$MOD_RISK_EDIT == "Adverse")], 
                                 "lower" = c1_mod$lower[which(c1_mod$MOD_RISK_EDIT == "Adverse")], 
                                 "upper" = c1_mod$upper[which(c1_mod$MOD_RISK_EDIT == "Adverse")]), 
          
          "Favorable" = data.frame("time" = c1_mod$time[which(c1_mod$MOD_RISK_EDIT == "Favorable")], 
                                   "censor" = c1_mod$n.censor[which(c1_mod$MOD_RISK_EDIT == "Favorable")], 
                                   "surv"  = c1_mod$surv[which(c1_mod$MOD_RISK_EDIT == "Favorable")], 
                                   "lower" = c1_mod$lower[which(c1_mod$MOD_RISK_EDIT == "Favorable")], 
                                   "upper" = c1_mod$upper[which(c1_mod$MOD_RISK_EDIT == "Favorable")]),
          
          "Intermediate" = data.frame("time" = c1_mod$time[which(c1_mod$MOD_RISK_EDIT == "Intermediate")], 
                                      "censor" = c1_mod$n.censor[which(c1_mod$MOD_RISK_EDIT == "Intermediate")], 
                                      "surv"  = c1_mod$surv[which(c1_mod$MOD_RISK_EDIT == "Intermediate")], 
                                      "lower" = c1_mod$lower[which(c1_mod$MOD_RISK_EDIT == "Intermediate")], 
                                      "upper" = c1_mod$upper[which(c1_mod$MOD_RISK_EDIT == "Intermediate")])
        )
        
      } else if(cen_type == "censct"){ # SCT as censored then take c3
        res_surv <- list(
          "Adverse" = data.frame("time" = c3_mod$time[which(c3_mod$MOD_RISK == "Adverse")], 
                                 "censor" = c3_mod$n.censor[which(c3_mod$MOD_RISK == "Adverse")], 
                                 "surv"  = c3_mod$surv[which(c3_mod$MOD_RISK == "Adverse")], 
                                 "lower" = c3_mod$lower[which(c3_mod$MOD_RISK == "Adverse")], 
                                 "upper" = c3_mod$upper[which(c3_mod$MOD_RISK == "Adverse")]), 
          
          "Favorable" = data.frame("time" = c3_mod$time[which(c3_mod$MOD_RISK == "Favorable")], 
                                   "censor" = c3_mod$n.censor[which(c3_mod$MOD_RISK == "Favorable")], 
                                   "surv"  = c3_mod$surv[which(c3_mod$MOD_RISK == "Favorable")], 
                                   "lower" = c3_mod$lower[which(c3_mod$MOD_RISK == "Favorable")], 
                                   "upper" = c3_mod$upper[which(c3_mod$MOD_RISK == "Favorable")]),
          
          "Intermediate" = data.frame("time" = c3_mod$time[which(c3_mod$MOD_RISK == "Intermediate")], 
                                      "censor" = c3_mod$n.censor[which(c3_mod$MOD_RISK == "Intermediate")], 
                                      "surv"  = c3_mod$surv[which(c3_mod$MOD_RISK == "Intermediate")], 
                                      "lower" = c3_mod$lower[which(c3_mod$MOD_RISK == "Intermediate")], 
                                      "upper" = c3_mod$upper[which(c3_mod$MOD_RISK == "Intermediate")])
        )
        
        res_surv_edit <- list(
          "Adverse" = data.frame("time" = c3_mod$time[which(c3_mod$MOD_RISK_EDIT == "Adverse")], 
                                 "censor" = c3_mod$n.censor[which(c3_mod$MOD_RISK_EDIT == "Adverse")], 
                                 "surv"  = c3_mod$surv[which(c3_mod$MOD_RISK_EDIT == "Adverse")], 
                                 "lower" = c3_mod$lower[which(c3_mod$MOD_RISK_EDIT == "Adverse")], 
                                 "upper" = c3_mod$upper[which(c3_mod$MOD_RISK_EDIT == "Adverse")]), 
          
          "Favorable" = data.frame("time" = c3_mod$time[which(c3_mod$MOD_RISK_EDIT == "Favorable")], 
                                   "censor" = c3_mod$n.censor[which(c3_mod$MOD_RISK_EDIT == "Favorable")], 
                                   "surv"  = c3_mod$surv[which(c3_mod$MOD_RISK_EDIT == "Favorable")], 
                                   "lower" = c3_mod$lower[which(c3_mod$MOD_RISK_EDIT == "Favorable")], 
                                   "upper" = c3_mod$upper[which(c3_mod$MOD_RISK_EDIT == "Favorable")]),
          
          "Intermediate" = data.frame("time" = c3_mod$time[which(c3_mod$MOD_RISK_EDIT == "Intermediate")], 
                                      "censor" = c3_mod$n.censor[which(c3_mod$MOD_RISK_EDIT == "Intermediate")], 
                                      "surv"  = c3_mod$surv[which(c3_mod$MOD_RISK_EDIT == "Intermediate")], 
                                      "lower" = c3_mod$lower[which(c3_mod$MOD_RISK_EDIT == "Intermediate")], 
                                      "upper" = c3_mod$upper[which(c3_mod$MOD_RISK_EDIT == "Intermediate")])
        )
      }
      
      t_max <- t_max 
      tux <- tux 
      
      if(jj == 1){
        t1 <- res_surv$Adverse$time[which(res_surv$Adverse$time <= t_max)]
        t2 <- res_surv$Favorable$time[which(res_surv$Favorable$time <= t_max)]
        t3 <- res_surv$Intermediate$time[which(res_surv$Intermediate$time <= t_max)]
        tu <- sort(unique(c(t1, t2, t3)))
        
        
        fitkm <- simple_survival_prediction_km(tu, tux, t_max, u1_surv = res_surv, u1 = KMdata, u1_test = KMdata_test, u1_sct = KMdata_sct, u1_sct_test = KMdata_sct_test, cen_type)
        t_risksetROC <- fitkm$t_risksetROC
        t_survivalROC <- fitkm$t_survivalROC
        ibrier <- fitkm$ibrier
        brier_365 <- unlist(fitkm$brier_365)
        
        ## Model specific prediction
        ## Model specific prediction
        fitrdg <- simple_survival_prediction(tu, var_model = c(vrm[vrm %!in% c("ELN_RISK_GROUP", "MOD_RISK_EDIT")]), u1 = KMdata, u1_test = KMdata_test, u1_sct = KMdata_sct, u1_sct_test = KMdata_sct_test, cen_type)
        t_risksetROC_ridge <- fitrdg$t_risksetROC_ridge
        t_survivalROC_ridge <- fitrdg$t_survivalROC_ridge
        ibrier_cridge <- fitrdg$ibrier_cridge
        median(t_survivalROC_ridge, na.rm = T)
      }
      
      if(jj == 2){
        t1_edit <- res_surv_edit$Adverse$time[which(res_surv_edit$Adverse$time <= t_max)]
        t2_edit <- res_surv_edit$Favorable$time[which(res_surv_edit$Favorable$time <= t_max)]
        t3_edit <- res_surv_edit$Intermediate$time[which(res_surv_edit$Intermediate$time <= t_max)]
        tu_edit <- sort(unique(c(t1_edit, t2_edit, t3_edit)))
        
        fitkm_edit <- simple_survival_prediction_km(tu = tu_edit, tux, t_max, u1_surv = res_surv_edit, u1 = KMdata %>% mutate(MOD_RISK = MOD_RISK_EDIT), u1_test = KMdata_test %>% mutate(MOD_RISK = MOD_RISK_EDIT), 
                                                    u1_sct = KMdata_sct %>% mutate(MOD_RISK = MOD_RISK_EDIT), u1_sct_test = KMdata_sct_test %>% mutate(MOD_RISK = MOD_RISK_EDIT), cen_type)
        t_risksetROC_edit <- fitkm_edit$t_risksetROC
        t_survivalROC_edit <- fitkm_edit$t_survivalROC
        ibrier_edit <- fitkm_edit$ibrier
        brier_365_edit <- unlist(fitkm_edit$brier_365)
        
        ## Here we directly use mod_risk_edit variable in the var_model and thus not needed to replace mod_risk with mod_risk_edit as above 
        fitrdg_edit <- simple_survival_prediction(tu = tu_edit, var_model = c(vrm[vrm %!in% c("ELN_RISK_GROUP", "MOD_RISK")]), u1 = KMdata, u1_test = KMdata_test, u1_sct = KMdata_sct, u1_sct_test = KMdata_sct_test, cen_type)
        t_risksetROC_ridge_edit <- fitrdg_edit$t_risksetROC_ridge
        t_survivalROC_ridge_edit <- fitrdg_edit$t_survivalROC_ridge
        ibrier_cridge_edit <- fitrdg_edit$ibrier_cridge
      }
      
      
      ## Model specific ELN prediction
      ## Model specific ELN prediction
      fitrdg_eln <- simple_survival_prediction(tu = tu, var_model = c(vrm[vrm %!in% c("MOD_RISK", "MOD_RISK_EDIT")], "ELN_RISK_GROUP"), u1 = KMdata, u1_test = KMdata_test, u1_sct = KMdata_sct, u1_sct_test = KMdata_sct_test, cen_type)
      t_risksetROC_ridge_eln <- fitrdg_eln$t_risksetROC_ridge
      t_survivalROC_ridge_eln <- fitrdg_eln$t_survivalROC_ridge
      ibrier_cridge_eln <- fitrdg_eln$ibrier_cridge
      
      ### KM specific prediction by ELN 
      if(cen_type == "all"){ # ITT
        res_surv_eln <- list(
          "Adverse" = data.frame("time" = c1_eln$time[which(c1_eln$ELN_RISK_GROUP == "Adverse")], 
                                 "censor" = c1_eln$n.censor[which(c1_eln$ELN_RISK_GROUP == "Adverse")],
                                 "surv"  = c1_eln$surv[which(c1_eln$ELN_RISK_GROUP == "Adverse")], 
                                 "lower" = c1_eln$lower[which(c1_eln$ELN_RISK_GROUP == "Adverse")], 
                                 "upper" = c1_eln$upper[which(c1_eln$ELN_RISK_GROUP == "Adverse")]), 
          
          "Favorable" = data.frame("time" = c1_eln$time[which(c1_eln$ELN_RISK_GROUP == "Favorable")], 
                                   "censor" = c1_eln$n.censor[which(c1_eln$ELN_RISK_GROUP == "Favorable")],
                                   "surv"  = c1_eln$surv[which(c1_eln$ELN_RISK_GROUP == "Favorable")], 
                                   "lower" = c1_eln$lower[which(c1_eln$ELN_RISK_GROUP == "Favorable")], 
                                   "upper" = c1_eln$upper[which(c1_eln$ELN_RISK_GROUP == "Favorable")]),
          
          "Intermediate" = data.frame("time" = c1_eln$time[which(c1_eln$ELN_RISK_GROUP == "Intermediate")], 
                                      "censor" = c1_eln$n.censor[which(c1_eln$ELN_RISK_GROUP == "Intermediate")],
                                      "surv"  = c1_eln$surv[which(c1_eln$ELN_RISK_GROUP == "Intermediate")], 
                                      "lower" = c1_eln$lower[which(c1_eln$ELN_RISK_GROUP == "Intermediate")], 
                                      "upper" = c1_eln$upper[which(c1_eln$ELN_RISK_GROUP == "Intermediate")])
        )
      } else if(cen_type == "censct"){ # SCT as censored then take c3
        res_surv_eln <- list(
          "Adverse" = data.frame("time" = c3_eln$time[which(c3_eln$ELN_RISK_GROUP == "Adverse")], 
                                 "censor" = c3_eln$n.censor[which(c3_eln$ELN_RISK_GROUP == "Adverse")],
                                 "surv"  = c3_eln$surv[which(c3_eln$ELN_RISK_GROUP == "Adverse")], 
                                 "lower" = c3_eln$lower[which(c3_eln$ELN_RISK_GROUP == "Adverse")], 
                                 "upper" = c3_eln$upper[which(c3_eln$ELN_RISK_GROUP == "Adverse")]), 
          
          "Favorable" = data.frame("time" = c3_eln$time[which(c3_eln$ELN_RISK_GROUP == "Favorable")], 
                                   "censor" = c3_eln$n.censor[which(c3_eln$ELN_RISK_GROUP == "Favorable")],
                                   "surv"  = c3_eln$surv[which(c3_eln$ELN_RISK_GROUP == "Favorable")], 
                                   "lower" = c3_eln$lower[which(c3_eln$ELN_RISK_GROUP == "Favorable")], 
                                   "upper" = c3_eln$upper[which(c3_eln$ELN_RISK_GROUP == "Favorable")]),
          
          "Intermediate" = data.frame("time" = c3_eln$time[which(c3_eln$ELN_RISK_GROUP == "Intermediate")], 
                                      "censor" = c3_eln$n.censor[which(c3_eln$ELN_RISK_GROUP == "Intermediate")],
                                      "surv"  = c3_eln$surv[which(c3_eln$ELN_RISK_GROUP == "Intermediate")], 
                                      "lower" = c3_eln$lower[which(c3_eln$ELN_RISK_GROUP == "Intermediate")], 
                                      "upper" = c3_eln$upper[which(c3_eln$ELN_RISK_GROUP == "Intermediate")])
        )
      }
      
      t_max <- t_max
      
      t1 <- res_surv_eln$Adverse$time[which(res_surv_eln$Adverse$time <= t_max)]
      t2 <- res_surv_eln$Favorable$time[which(res_surv_eln$Favorable$time <= t_max)]
      t3 <- res_surv_eln$Intermediate$time[which(res_surv_eln$Intermediate$time <= t_max)]
      tu_eln <- sort(unique(c(t1, t2, t3)))
      tux_eln <- tux
      
      fitkm_eln <- simple_survival_prediction_eln_km(tu = tu_eln, tux_eln, t_max, u1_surv =  res_surv_eln, u1 = KMdata, u1_test = KMdata_test, u1_sct = KMdata_sct, u1_sct_test = KMdata_sct_test, cen_type) 
      t_risksetROC_eln <- fitkm_eln$t_risksetROC
      t_survivalROC_eln <- fitkm_eln$t_survivalROC
      ibrier_eln <- fitkm_eln$ibrier
      brier_365_eln <- unlist(fitkm_eln$brier_365)
    } 
  } ## loop for jj ends 
      
      
      # Create a list of vectors with varying lengths
      list_of_vectors <- list("Strongly adverse" = sa, "Moderately adverse" = ma, "Moderately favorable" = mf, "Strongly favorable" = sf) # method-1
      list_of_vectors_edit <- list("Strongly adverse" = sa_edit, "Moderately adverse" = ma_edit, "Moderately favorable" = mf_edit, "Strongly favorable" = sf_edit) # method-2
      
      # Find the maximum length among the vectors
      max_length <- max(lengths(list_of_vectors))
      max_length_edit <- max(lengths(list_of_vectors_edit))
      
      # Pad the vectors with NA to make them of equal length
      padded_list <- lapply(list_of_vectors, function(x) {
        length(x) <- max_length
        x
      })
      padded_list_edit <- lapply(list_of_vectors_edit, function(x) {
        length(x) <- max_length_edit
        x
      })
      
      # Create a data frame from the padded list
      risk_df <- as.data.frame(padded_list)
      risk_df_edit <- as.data.frame(padded_list_edit)
      
      test_risk <- KMdata_test$MOD_RISK
      test_risk_edit <- KMdata_test$MOD_RISK_EDIT
      
      test_mrn <- KMdata_test$MRN
      
      
      saveRDS(
      list(sa = sa, ma = ma, na = na, sf = sf, mf = mf, nf = nf, risk_df = risk_df, test_risk = test_risk, test_mrn = test_mrn, tu = tu, tux = tux,  
           res_surv_eln = res_surv_eln,  
           res_surv = res_surv , t_risksetROC = t_risksetROC , t_survivalROC = t_survivalROC , ibrier = ibrier , brier_365 = brier_365, 
           t_risksetROC_eln = t_risksetROC_eln, t_survivalROC_eln = t_survivalROC_eln, ibrier_eln = ibrier_eln, brier_365_eln = brier_365_eln,
           
           ibrier_cridge = ibrier_cridge, t_survivalROC_ridge = t_survivalROC_ridge, t_risksetROC_ridge = t_risksetROC_ridge, 
           ibrier_cridge_eln = ibrier_cridge_eln, t_survivalROC_ridge_eln = t_survivalROC_ridge_eln, t_risksetROC_ridge_eln = t_risksetROC_ridge_eln,
           
           sa_edit = sa_edit, ma_edit = ma_edit, na_edit = na_edit, sf_edit = sf_edit, mf_edit = mf_edit, nf_edit = nf_edit, risk_df_edit = risk_df_edit, test_risk_edit = test_risk_edit, 
           res_surv_edit = res_surv_edit, t_risksetROC_edit = t_risksetROC_edit, t_survivalROC_edit = t_survivalROC_edit, ibrier_edit = ibrier_edit, 
           brier_365_edit = brier_365_edit, ibrier_cridge_edit = ibrier_cridge_edit, t_survivalROC_ridge_edit = t_survivalROC_ridge_edit, t_risksetROC_ridge_edit = t_risksetROC_ridge_edit), 
      file = paste0("RISK_DF_FC_", fc_include,"_CD71_" , exclude_cd71, "_", mis_type, "_", cen_type,"_", good_def,".rds"))
}




## survival test 
# lsurv = survival model fit 
# test1 = Type of test (e.g., equality) 
# test2 = Method of tests (e.g., LR = Log-rank)
# test_p = Weights for weighted LR tests (FH = Fleming harrington)
# test_q = Weights for weighted LR tests (FH = Fleming harrington)
# inx = Index for test-statistics
survtest <- function(lsurv, test1, test2, test_p, test_q, inx){
  
  survMisc::comp(lsurv, p = test_p[1], q = test_p[1], reCalc = TRUE)
  
  if(test1 == "Equality"){
    pval <- attr(lsurv, "lrt")
    chi <- 1 - pchisq(pval$chiSq, pval$df)[inx]  
  }
  
  if(test1 == "Ordered differences"){
    pval <- attr(lsurv, "tft")
    chi <- 1 - pchisq(pval$tft$chiSq, pval$tft$df)[inx]
  }
  
 
  result  <- data.frame("Test" =  ifelse(test2 == "Log-rank", "LR",
                                         ifelse(test2 == "Gehan-Breslow", "GB",
                                                ifelse(test2 == "Tarone-Ware", "TW",
                                                       ifelse(test2 == "Peto-Prentice", "PP",
                                                              ifelse(test2 == "Peto-Prentice-Anderson", "PPA",
                                                                     ifelse(test2 == "Fleming-Harrington", "FH", "LR"))))))
                        , "p-value" = chi)
  
  
  if("Fleming-Harrington" %in% test2){
    survMisc::comp(lsurv, p = test_p[which(test2 == "Fleming-Harrington")],
                   q = test_q[which(test2 == "Fleming-Harrington")], reCalc = TRUE)
    
    if(test1 == "Equality"){
      pval <- attr(lsurv, "lrt")
      chi <- 1 - pchisq(pval$chiSq, pval$df)[6]  
    }
    
    if(test1 == "Ordered differences"){
      pval <- attr(lsurv, "tft")
      chi <- 1 - pchisq(pval$tft$chiSq, pval$tft$df)[6]
    }
    
    result$p.value[which(result$Test == "FH")] <- chi
  }
  result
}

survtest_renyi <- function(lsurv, test1, test2, test_p, test_q, inx){
  comp(lsurv, p = 0, q = 0)
  if(test1 == "Equality"){
    pval <- attr(lsurv, "sup")
    chi <- pval$pSupBr[inx]  
  }
  
  result  <- data.frame("Test" = test2, "p-value" = chi)
  result
}


## Variable adjusted survival analysis
#
# conf_int = Boolean whether 95% CI needs to be generated
# median_line = Boolean whether median survival days needs to be reported
# col_set = Color for curves
# name_var = Variable name how it needs to be displayed in the plot
# xsurv = original variable name in dataset
# Main = Analytical dataset; both time and censor needs to be precomputed and added in this dataset
# font.ytickslab, font.xtickslab, font.title, font.x, font.y, censor.size, xlab, conf.int.alpha, 
# conf.int.style, legend_name, wghts are all graphical paramaters 
plotInput2_paperML <- function(conf_int = TRUE, median_line = TRUE,
                               test1 = "Equality", test2 = "Log-rank",
                               test_p = 1, # weight for Grehan-breslow
                               test_q = 1, # weight for Grehan-breslow 
                               col_set,    # set of colors
                               name_var = "CD7",
                               xsurv = "FC_CD_OTHER", 
                               xsurv2 = "none",
                               xsurv3 = "none",
                               Main, 
                               font.ytickslab = 17, font.xtickslab = 17,
                               font.title = 1,
                               font.x = 15, font.y  = 15, censor.size = 5,
                               xlab = "Time (days) post-treatment",
                               conf.int.alpha = 0.2, conf.int.style = "ribbon",
                               legend_name = NULL,
                               wghts = NULL){
  
  if(xsurv != "none" | xsurv2 != "none" | xsurv3 != "none"){
    
    hh <- c(xsurv, xsurv2, xsurv3)
    hh_del <- list()
    
    aa <- 1.5
    dplot2 <- Main
    
    if(nrow(dplot2) == 0){ 
      stop(paste0("We can not perform analysis as all values for this variable are Not performed!"))
    }
    
    if(xsurv != "none" & xsurv2 == "none" & xsurv3 == "none"){
      
      ## deleting all unknowns 
      dplot2 <- subdataset <- dplot2[, colnames(dplot2) %in% c(xsurv, "time", "censor", "wghts")] %>% #mutate_if(is.character, as.factor) %>%
        mutate(across(
          where(is.factor),
          ~if_else(. %in% c("Missing", "UNKNOWN", "NP", "Not performed"), NA, .)
        )
        ) %>% drop_na()
      
      
      AA <- xsurv
      Category <- paste0("droplevels(as.factor(dplot2$", AA, "))", sep = "")
      Category <- eval(parse(text = Category))
      table(Category)
      
      if( any(c("Y", "N") %in% unique(Category, na.rm = T)) == TRUE ){ # if atleast one "Y" or "N" then do following
        Category <- gsub("Y", "Yes", Category) # replacing if Y is there
        Category <- gsub("N", "No", Category)  # replacng if N is there
      }
      
      ut <- unique(sort(Category, na.rm = NA))
      lt <- length(unique(Category, na.rm = NA))
      
      
      
      if(conf_int ==  TRUE){
        conf = TRUE
      } else{
        conf = FALSE
      }
      
      if(median_line == TRUE){
        median = TRUE
      } else{
        median = FALSE
      }
      
      prnt <- plotFunc1_paperML(dplot2sub = dplot2, AA, Category, ut, lt, col_set, conf, median,
                              test1 = test1, 
                              test2 = test2, 
                              test_p = as.numeric(test_p),
                              test_q = as.numeric(test_q),
                              name_var,
                              font.ytickslab = font.ytickslab, font.xtickslab = font.xtickslab,
                              font.title = font.title,
                              font.x = font.x, font.y  = font.y, censor.size = censor.size,
                              xlab = xlab,
                              conf.int.alpha = conf.int.alpha, conf.int.style = conf.int.style,
                              legend_name = legend_name, wghts = wghts)
    }
    
   
  } 
  output = list(prnt = prnt, subdataset = subdataset)
  
}  


## Adjusted KM (one covar)
plotFunc1_paperML <- function(dplot2sub, AA, Category, ut, lt, col_set, conf, median,
                              test1, test2, test_p, test_q, name_var,
                              font.ytickslab = 17, font.xtickslab = 17,
                              font.title = 1,
                              font.x = 15, font.y  = 15, censor.size = censor.size,
                              xlab = "Time (days) post-treatment",
                              conf.int.alpha = 0.2, conf.int.style = "ribbon", 
                              legend_name = NULL, wghts = NULL){
  
  aa <- 1.1; cexcol <- 0.7; cexrow <- 0.6; base_tt <- 12; base_gg <- 12; cex_gg <- 0.9
  
  Status <- ifelse(dplot2sub$censor == 1, "Deceased", "Alive")
    
  ## Dataframe creation
  QQ <- paste0("", AA, " <- Category" )
  QQ <- eval(parse(text = QQ))
  
  
  tt <-  paste0("as.data.frame(table(", AA, ", Status))", sep = "")
  tt <- eval(parse(text = tt))
  
  pos <- unlist(lapply(seq_len(lt), function(vv) c(vv, (vv + lt))))
  
  if(length(unique(tt$Status)) != 1){ # There are both deceased & alive
    tt$Size <- rep(unlist(lapply(seq_len(lt), function(vv) sum(tt$Freq[c(vv, (vv + lt))]))), 2)
  } else{ # There is only alive
    tt$Size <- tt$Freq
  }
  
  QQ <- paste0("tt <- tt[c('", AA, "', 'Size', 'Status',  'Freq')]")
  QQ <- eval(parse(text = QQ))
  
  colnames(tt)[(ncol(tt) - 2)] <- "Total N"
  colnames(tt)[ncol(tt)] <- "Frequency"
  
  
  if(nrow(tt) == 1 & tt$Status[1] == "Deceased"){
    stop(paste0("KM can not be estimated! There is only one subject (i.e., MRN : ", dplot2sub$MRN, ") & who is deceased!"))
  }
  
  fc <- colorRampPalette(c("gray100", "gray65"))
  fc <- fc(24) # creating 16 sheds
  
  
  fills <- matrix(rep(fc[1 : lt], 2 * ncol(tt)), nrow = 2 * lt, byrow = F)  #  2 =  dead/alive & 5 columns: Key/cat1/stat/freq/Size
  cols <- matrix(rep(rep("black", lt), 2 * ncol(tt)), nrow = 2 * lt, byrow = F)     # color matrix for text
  
  tt <- tt[pos, ]; fills <- fills[pos, ]; cols <- cols[pos, ] # merging alive+deceased together
  cols[seq(2, 2 * lt, 2), (ncol(tt) - 2)] <-  fills[seq(2, 2 * lt, 2), (ncol(tt) - 2)] # replacing black text with row color to hide
  
  if(length(which(is.na(tt$`Total N`))) != 0){ # deleting extra row that are created at previous step | MUST BE NUMERIC/ FOR CLASS IT TAKES NA AS ANOTHER CLASS/LABEL
    extra_row <- which(is.na(tt$Frequency)) 
    tt <- tt[- extra_row, ]
    cols <- cols[- extra_row, ]
    fills <- fills[- extra_row, ]
  }
  
  
  rownames(tt) <- seq_len(nrow(tt))
  
  g <- tableGrob(tt, theme = ttheme_minimal(
    core = list(bg_params = list(fill = fills, col = NA),
                fg_params = list(fontface = 3, col = cols)), base_size = base_tt,
    colhead = list(fg_params = list(rot = 55, col = "black", fontface = 4L, cex = cexcol)),
    rowhead = list(fg_params = list(col = "black", fontface = 3L, cex = cexrow))
  ))
  
  
  separators <- replicate(ncol(g) - 2, segmentsGrob(x1 = unit(0, "npc"), gp = gpar(lty = 2)), simplify = FALSE)
  g <- gtable::gtable_add_grob(g, grobs = separators, t = 2, b = nrow(g), l = seq_len(ncol(g) - 2) + 2)
  
  
  
  ## KM curves
  if(length(wghts) == 0){ 
    wghts <- NULL
    lsurv2 <- paste0("survfit(Surv(time, censor, type = 'right') ~ ", AA,  ", data = dplot2sub, type = 'kaplan-meier', error = 'tsiatis') ")
    lsurv2 <- eval(parse(text = lsurv2))
  } else{
    wghts <- dplot2sub$weights
    lsurv2 <- paste0("survfit(Surv(time, censor, type = 'right') ~ ", AA,  ", data = dplot2sub, type = 'kaplan-meier', error = 'tsiatis', weights = wghts) ")
    lsurv2 <- eval(parse(text = lsurv2))
  }
  
  lt <- paste0("length(unique(dplot2sub$", AA, "))", sep = "")
  lt <- eval(parse(text = lt))
  
  
  brk <- round(seq(1, (max(dplot2sub$time)- 60), length.out = 6)[2] - seq(1, (max(dplot2sub$time)- 60), length.out = 6)[1], 0)
  if(brk < 0){
    brk <- round(seq(1, (max(dplot2sub$time)- 1), length.out = 6)[2] - seq(1, (max(dplot2sub$time)- 1), length.out = 6)[1], 0)  
  }
  
  
  if(lt != 1 & length(unique(tt$Status)) != 1  | # length(unique(tt$`Total N`)) != 1
     
     length(which(ut %in% c("Experienced CR/CRi", "Experienced CR", "Other", "Yes", "No") == TRUE)) > 1){ # Test can not be performed if we have single category
    
    if(length(wghts) == 0){ 
      lsurv <- paste0("ten(Surv(time, censor, type = 'right') ~ ", AA,  ", data = dplot2sub, type = 'kaplan-meier', error = 'tsiatis') ")
      lsurv <- eval(parse(text = lsurv))
    } else{
      lsurv <- paste0("ten(Surv(time, censor, type = 'right') ~ ", AA,  ", data = dplot2sub, type = 'kaplan-meier', error = 'tsiatis', weights = wghts) ")
      lsurv <- eval(parse(text = lsurv))
    }
    
    inx <- ifelse(test2 == "Log-rank", 1,
                  ifelse(test2 == "Gehan-Breslow", 2,
                         ifelse(test2 == "Tarone-Ware", 3,
                                ifelse(test2 == "Peto-Prentice", 4,
                                       ifelse(test2 == "Peto-Prentice-Anderson", 5,
                                              ifelse(test2 == "Fleming-Harrington", 6, 1))))))
    
  }
  
  med_val <- NULL 
  if(median == FALSE){
    if(conf == FALSE){
      
      g1 <- ggsurvplot(lsurv2, title = "",
                       linetype = "solid", legend = "none", 
                       palette = col_set[1 : lt], 
                       font.ytickslab = c(font.ytickslab, "plain", "black"),
                       font.xtickslab = c(font.xtickslab, "plain", "black"),
                       censor.size = censor.size, tables.y.text = FALSE,
                       ggtheme = theme_bw(),
                       xlab = xlab,
                       font.title = c(font.title, "bold", "blue"),
                       font.x = c(font.x, "plain", "black"),
                       font.y = c(font.y, "plain", "black"),
                       tables.theme = clean_theme(), data = dplot2sub,
                       risk.table = "nrisk_cumcensor")
    } else{
      
      g1 <- ggsurvplot(lsurv2, title = "",
                       conf.int = TRUE, conf.int.style = conf.int.style , conf.int.alpha = conf.int.alpha ,
                       linetype = "solid", legend = "none",
                       palette = col_set[1 : lt], 
                       font.ytickslab = c(font.ytickslab, "plain", "black"),
                       font.xtickslab = c(font.xtickslab, "plain", "black"),
                       censor.size = censor.size, tables.y.text = FALSE,
                       ggtheme = theme_bw(),
                       xlab = xlab,
                       font.title = c(font.title, "bold", "blue"),
                       font.x = c(font.x, "plain", "black"),
                       font.y = c(font.y, "plain", "black"),
                       tables.theme = clean_theme(), data = dplot2sub,
                       risk.table = "nrisk_cumcensor")
    }
  }
  
  
  if(median == TRUE){
    if(conf == FALSE){
      g1 <- ggsurvplot(lsurv2, title = "",
                       surv.median.line = "hv",
                       linetype = "solid", legend = "none",
                       palette = col_set[1 : lt], 
                       font.ytickslab = c(font.ytickslab, "plain", "black"),
                       font.xtickslab = c(font.xtickslab, "plain", "black"),
                       censor.size = censor.size, tables.y.text = FALSE,
                       ggtheme = theme_bw(),
                       xlab = xlab,
                       font.title = c(font.title, "bold", "blue"),
                       font.x = c(font.x, "plain", "black"),
                       font.y = c(font.y, "black", "black"),
                       tables.theme = clean_theme(), data = dplot2sub,
                       risk.table = "nrisk_cumcensor")
      g1$plot <- g1$plot
    } else{
      g1 <- ggsurvplot(lsurv2, title = "",
                       surv.median.line = "hv",
                       conf.int = TRUE, conf.int.style = conf.int.style , conf.int.alpha = conf.int.alpha ,
                       linetype = "solid", legend = "none",
                       palette = col_set[1 : lt], 
                       font.ytickslab = c(font.ytickslab, "plain", "black"),
                       font.xtickslab = c(font.xtickslab, "plain", "black"),
                       censor.size = censor.size, tables.y.text = FALSE,
                       ggtheme = theme_bw(),
                       xlab = xlab,
                       font.title = c(font.title, "plain", "black"),
                       font.x = c(font.x, "plain", "black"),
                       font.y = c(font.y, "plain", "black"),
                       tables.theme = clean_theme(), data = dplot2sub,
                       risk.table = "nrisk_cumcensor")
    }
    
    
    # Display median
    med_val <- as.numeric(surv_median(lsurv2)[, 2])
    for(jj in 1 : length(med_val)){
      
      g1$plot <- g1$plot + ggplot2::annotate("text",
                                             x = med_val[jj], y = 0,   col = "red", angle = 90, fontface = 2,  
                                             label = paste("", med_val[jj], ""), size = 5)
    }
    
  } # MEDIAN true END
  
  g5 <- g1$plot
  g6 <- pv <- pframe <- NULL
  
  if(test1 == "Equality" | test1 == "Ordered differences"){
    
    if(lt != 1 & length(unique(tt$Status)) != 1 | 
       
       length(which(ut %in% c("Experienced CR/CRi", "Experienced CR", "Other", "Yes", "No") == TRUE)) > 1) {
      iu <- paste0("which(table(dplot2sub$", AA, ") == 1)") # Need to have at-least 2 cases for testing | therefore dropping that group
      iu <- eval(parse(text = iu))
      
      if(length(iu) > 0){
        # dropping the class that has only 1 count in testing
        iu2 <- paste0("which(dplot2sub$", AA, "== names(table(dplot2sub$", AA,")[iu]))")
        iu2 <- eval(parse(text = iu2))
        dplot2sub2 <- dplot2sub[- iu2, ]
      } else{
        dplot2sub2 <- dplot2sub
      }
      if(length(wghts) == 0){ 
        lsurv <- paste0("ten(Surv(time, censor, type = 'right') ~ ", AA,  ", data = dplot2sub, type = 'kaplan-meier', error = 'tsiatis') ")
        lsurv <- eval(parse(text = lsurv))
      } else{
        lsurv <- paste0("ten(Surv(time, censor, type = 'right') ~ ", AA,  ", data = dplot2sub, type = 'kaplan-meier', error = 'tsiatis', weights = wghts) ")
        lsurv <- eval(parse(text = lsurv))
      }
      
      testres <- survtest(lsurv, test1, test2, test_p, test_q, inx = inx)
      testres_fh <- survtest(lsurv, test1, test2, test_p = c(0, 0, 0), test_q = c(0, 0, 1), inx = inx) # always fleming-harrington
      testres <- rbind(testres, testres_fh[which(testres_fh[, 1] %in% c("FH", "Fleming-Harrington")), ])
      
      re <- which(testres$Test %in% "FH")
      testres$Test <- as.character(testres$Test)
      
      if(length(re) == 1){
        testres$Test[re[1]]  <- "FH:e"
      }
      if(length(re) == 2){
        testres$Test[re[1]]  <- "FH:e"
        testres$Test[re[2]]  <- "FH:l"
      }
      re <- which(testres$Test %in% "LR")
      if(length(re) > 1){
        testres <- testres[-re[-1], ]
      }
      
      
      g6 <- tableGrob(testres, rows = NULL, theme = ttheme_minimal(
        core = list(fg_params = list(fontface = 3, col = "firebrick2")), base_size = base_gg,
        colhead = list(fg_params = list(col = "black", fontface = 4L, cex = cex_gg))))
      
      g6 <- g5 + annotation_custom(grob = g6, xmin = 0, ymin = 0, ymax = 0.15)
      
      pv <- ifelse(sprintf("%0.4f", testres$p.value) == "0.0000", "<0.0001", sprintf("%0.4f", testres$p.value))
      
      e1 <- c(testres$Test); e2 <- pv
      
      # mAX-combo
      if(test1 == "Equality" & lt == 2 & "Max-Combo" %in% test2){
        
        grp <- paste0("dplot2sub2$", AA, "")
        grp <- eval(parse(text = grp))
        set.seed(1212)
        try(mxc <- nph::logrank.maxtest(
          time = dplot2sub2$time,
          event = dplot2sub2$censor,
          group = grp, #
          alternative = c("two.sided"),
          rho = c(0, 0, 1),
          gamma = c(0, 1, 0),
          event_time_weights = NULL,
          algorithm = mvtnorm::GenzBretz(maxpts = 50000, abseps = 1e-05, releps = 0)
        ))
        
        if(class(mxc) %!in% "try-error"){
          pv <- c(pv, ifelse(sprintf("%0.4f", mxc$p.Bonf) == "0.0000", "<0.0001", sprintf("%0.4f", mxc$p.Bonf)))
          
          g1$plot <- g5 + ggplot2::annotate("text", x = -Inf, y = 0.80, hjust = -3.8, vjust = 2,  
                                            col = "black", label = paste0("MaxC   ", pv[length(pv)], ""), size = 5) 
          e1 <- c(e1, "MaxC"); e2 <- c(e2, pv[length(pv)])
        }
      }
      
      
      # Renyi
      if(test1 == "Equality" & lt == 2  & "Renyi" %in% test2){ 
        
        ## Renyi
        set.seed(1212)
        try(testres_renyi <- survtest_renyi(lsurv, test1, test2, test_p, test_q, inx = inx))
        if(class(testres_renyi) %!in% "try-error"){
          testres <- rbind(testres, testres_renyi[which(as.character(testres_renyi[, 1]) %in% c("LR", "Log-rank")), ])
          testres$Test[which(as.character(testres$Test) %in% c("LR", "Log-rank"))] <- "LR"
          testres$Test[tail(which(as.character(testres$Test) %in% c("LR", "Log-rank")), 1)] <- "LR:r" 
          pv <- ifelse(sprintf("%0.4f", testres$p.value) == "0.0000", "<0.0001", sprintf("%0.4f", testres$p.value))
          
          g1$plot <- g1$plot +  ggplot2::annotate("text", x = -Inf, y = 0.85, hjust = -4.4, vjust = 2, 
                                                  col = "black", label = paste0("", testres$Test[5],"  ", pv[length(pv)], ""), size = 5) # last one is added as Renyi
          
          e1 <- c(e1, testres$Test[5]); e2 <- c(e2, pv[length(pv)])
        }
      }
      
      ## mdir
      if(test1 == "Equality" & lt == 2  & "mdir2" %in% test2){ 
        drt <- data.frame("time" = dplot2sub$time, "event" = dplot2sub$censor, "group" = grp)
        ## mdir.tes
        set.seed(1212)
        try(system.time(try(a <- mdir.logrank(data = drt, rg = list(c(2, 2)), nperm = 10^3))))
        if(class(a) %!in% "try-error"){
          pv <- c(pv, a$p_value$Perm)
          a3 <- ifelse(sprintf("%0.4f", a$p_value$Perm) == "0.0000", "<0.0001", sprintf("%0.4f", a$p_value$Perm))
          g1$plot <- g1$plot +  ggplot2::annotate("text", x = -Inf, y = 0.75, hjust = -3.8, vjust = 2, 
                                                  col = "black", 
                                                  label = paste0("mDir2   ", a3, ""), size = 5)
          
          e1 <- c(e1, "mDir2"); e2 <- c(e2, a3)
        }
        
      }      
      
      ## mdir3
      if(test1 == "Equality" & lt == 2  & "mdir4" %in% test2){ 
        drt <- data.frame("time" = dplot2sub$time, "event" = dplot2sub$censor, "group" = grp) 
        ## mdir.tes
        set.seed(1212)
        try(system.time(try(a <- mdir.logrank(data = drt, rg = list(c(4, 4)), nperm = 10^3))))
        if(class(a) %!in% "try-error"){
          pv <- c(pv, a$p_value$Perm)
          a3 <- ifelse(sprintf("%0.4f", a$p_value$Perm) == "0.0000", "<0.0001", sprintf("%0.4f", a$p_value$Perm))
          g1$plot <- g1$plot +  ggplot2::annotate("text", x = -Inf, y = 0.70, hjust = -3.8, vjust = 2, 
                                                  col = "black", 
                                                  label = paste0("mDir4   ", a3, ""), size = 5)
          
          e1 <- c(e1, "mDir4"); e2 <- c(e2, a3)
        }
      }      
      
      ## rmst2
      if(test1 == "Equality" & lt == 2  & "RMST" %in% test2){ 
        
        a = rmst2(time = dplot2sub$time, status = dplot2sub$censor, arm = ifelse(as.character(grp) == ut[1], 0, 1), tau = NULL)
        a1 <- a$unadjusted.result[1, ]
        tau <- a$tau 
        
        a2 <- a1[names(a1) %in% "p"]
        pv <- c(pv, a2)
        a3 <- ifelse(sprintf("%0.4f", a2) == "0.0000", "<0.0001", sprintf("%0.4f", a2))
        
        g1$plot <- g1$plot +  ggplot2::annotate("text", x = -Inf, y = 0.65, hjust = -3.9, vjust = 2, 
                                                col = "black", label = paste0("RMST  ", a3, ""), size = 5) # last one is added as Renyi
        e1 <- c(e1, "RMST"); e2 <- c(e2, a3)
        
      }          
      
      ## KONP
      if(test1 == "Equality" & lt == 2  & "KONP" %in% test2){ 
        
        try(a <- try(konp_test(time = dplot2sub$time, 
                               status = dplot2sub$censor, 
                               group = ifelse(as.character(grp) == ut[1], 0, 1), 
                               n_perm = 10^3)))
        if(class(a) %!in% "try-error"){
          pv <- c(pv,  a$pv_chisq)
          a3 <- ifelse(sprintf("%0.4f", a$pv_chisq) == "0.0000", "<0.0001", sprintf("%0.4f", a$pv_chisq))
          g1$plot <- g1$plot +  ggplot2::annotate("text", x = -Inf, y = 0.60, hjust = -4, vjust = 2, 
                                                  col = "black", 
                                                  label = paste0("KONP  ", a3, ""), size = 5)
          e1 <- c(e1, "KONP"); e2 <- c(e2, a3)
        }
      }      
      
      pframe <- data.frame("name" = e1, "pval" = e2)
      
      
      g1$plot <- g1$plot + ggplot2::annotate("text", x = -Inf, y = Inf, hjust = -5.05, vjust = 2, 
                                             col = "black", label = paste0("", testres$Test[1],"  ", pv[1], ""), size = 5) +
        ggplot2::annotate("text", x = -Inf, y = 1, hjust = -4.8, vjust = 2,  
                          col = "black", label = paste0("", testres$Test[2],"  ", pv[2], ""), size = 5) +
        ggplot2::annotate("text", x = -Inf, y = 0.95, hjust = -4.3, vjust = 2, 
                          col = "black", label = paste0("", testres$Test[3],"  ", pv[3], ""), size = 5) +
        ggplot2::annotate("text", x = -Inf, y = 0.90, hjust = -4.5, vjust = 2, 
                          col = "black", label = paste0("", testres$Test[4],"  ", pv[4], ""), size = 5) 
      
      g1$plot 
      
    }
  } 
  
  g1$table <- g1$table + theme(axis.text.x = element_text(face = "plain", color = "black",
                                                          size = font.xtickslab),
                               plot.title = element_text(size = font.y, face = "bold", color = "black"))
  
  g3 <- as.ggplot(function(){
    plot(1, type = "n", xlab = "", ylab = "", ylim = c(0, 1), xlim = c(0, 2000), axes = FALSE)
    if(is.null(legend_name) == TRUE){
      legend(0, 1, title = "", 
             c(name_var, unlist(lapply(seq_len(lt), function(ff){
               vv <- paste0(" ' ", ut[ff], " ' ")
               vv <- eval(parse(text = vv))
             }))), cex = aa, horiz = FALSE,
             lwd = 7, lty = 1, bty = "n", col = c(NA, col_set[1 : lt]))
    } else {
      legend(0, 1, title = "", 
             c(name_var, legend_name), cex = aa, horiz = FALSE,
             lwd = 7, lty = 1, bty = "n", col = c(NA, col_set[1 : lt]))
    }
  })
  
  output = list(g5 = g5, g3 = g3, g1 = g1, g = g, g6 = g6, pv = pv, med_val = med_val, ut = ut, pframe = pframe)
}



## Predicting OOS by simple (only few variables) survival prediction function based on Cox-Ridge
# tu = Time at which evaluation needs to be performed
# var_model = Covariates to be adjusted in mCOXr
# u1 = Dataframe for full analytical set
# u1_test = Dataframe for test set
# u1_sct = Dataframe for full analytical set treating allo-HCT as censored
# u1_sct_test = Dataframe for test set treating allo-HCT as censored
# cen_type = Methods to deal with allo-HCT patients
simple_survival_prediction <- function(tu, var_model, u1, u1_test, u1_sct, u1_sct_test, cen_type) {
  
  if(cen_type ==  "all"){
    data <- u1 %>% filter(!is.na(censor))
    data_test <- u1_test %>% filter(!is.na(censor))
  } else if(cen_type == "censct"){
    data <- u1_sct %>% filter(!is.na(censor))
    data_test <- u1_sct_test %>% filter(!is.na(censor))
  }
  
  train_data <- data[, colnames(data) %in% c("time", "censor", var_model)]
  test_data <- data_test[, colnames(data_test) %in% c("time", "censor", var_model)]
  
  GG <- dplyr::bind_rows(test_data, train_data)
  nm0 <- paste0("", var_model, "", collapse = " + " )
  part1 <- paste0("Surv(time, censor) ~ ", nm0, "", collapse = " + ")
  
  fom <- paste0("survival::coxph(", part1, ", ties = 'breslow', iter.max = 1000, 
                  outer.max = 500, robust = FALSE,
                  singular.ok = TRUE, eps = 1e-1,
                  toler.chol = .Machine$double.eps^.15,
                  data = GG)")
  try(fom <- eval(parse(text =  fom)))
  ix_test <- seq_len(nrow(test_data))
  
  # extract model matrix
  mf <- model.matrix(fom)
  mf_train <- model.matrix(fom)[-ix_test, ]; dim(mf_train)
  mf_test <- as.matrix(model.matrix(fom)[ix_test, ]); dim(mf_test)  
  time <- train_data$time; censor <- train_data$censor
  yss <- survival::Surv(time, censor)
  
  centered = FALSE 
  set.seed(3000) 
  system.time(cv.fit <- cv.glmnet(mf_train, yss, family = "cox", nfolds = 10, type.measure = "deviance", alpha = 0,   
                                  gamma = seq(0, 1, length.out = 50), relax = FALSE, parallel = FALSE))
  lp <- as.numeric(predict(cv.fit, newx = data.matrix(mf_train), s = cv.fit$lambda.min, type = "link"))
  
  time <- train_data$time
  censor <- train_data$censor 
  t.unique <- tu
  
  alpha <- length(t.unique)
  for (i in 1L:length(t.unique)) {
    alpha[i] <- sum(time[censor == 1L] == t.unique[i]) / sum(exp(lp[time >= t.unique[i]]))
  }
  
  obj <- approx(t.unique, cumsum(alpha), yleft = 0, xout = t.unique, rule = 2)
  
  if (centered) {
    obj$y <- obj$y * exp(mean(lp))
  }
  obj$z <- exp(-obj$y)
  names(obj) <- c("times", "cumulative_base_hazard", "base_surv")
  
  
  lp <- as.numeric(predict(cv.fit, newx = data.matrix(mf), s = cv.fit$lambda.min, type = "link")) # with test data
  p <- exp(exp(lp) %*% (-t(obj$cumulative_base_hazard)))
  
  opt_param <- cv.fit$lambda.min
  stack_surv <- p[ix_test, ]     # predicted probabilities of all (NOTE: It gives predicted probabilities for both test and train data and first rows are test data)
  predicted_time <- lp[ix_test]
  
  Cstat_ridge <- rcorr.cens(-predicted_time, Surv(data_test$time, data_test$censor))
  
  ibrier_cridge <- SurvMetrics::IBS(Surv(data_test$time, data_test$censor), sp_matrix = stack_surv, t.unique)
  
  
  ## Cumulative case dynamic ROC
  t_survivalROC_ridge <- unlist(lapply(seq_len(length(t.unique)), function(tt) {
    survivalROC(Stime = data_test$time, 
                status = data_test$censor,
                marker = predicted_time,
                predict.time = t.unique[tt],
                method = "KM", # nearest neighbour estimation (NNE) & alternative is KM
                span = 0.25 * nrow(data_test)^(-0.20))$AUC # span don't need for KM
  }))
  
  t_risksetROC_ridge <- unlist(lapply(seq_len(length(t.unique)), function(tt) {
    risksetROC::risksetROC(Stime = data_test$time, 
                           status = data_test$censor,
                           marker = predicted_time,
                           predict.time = t.unique[tt], plot = F,
                           method = "Cox", # schoenfeld alternatives are Cox-PH
                           span = 0.25 * nrow(data_test)^(-0.20))$AUC # span don't need for schoenfeld
  }))
 
  output = list(ibrier_cridge = ibrier_cridge, t_survivalROC_ridge = t_survivalROC_ridge, t_risksetROC_ridge = t_risksetROC_ridge, Cstat_ridge = Cstat_ridge)
}


# simple survival prediction function based on KM 
simple_survival_prediction_km <- function(tu, tux, t_max, u1_surv, u1, u1_test, u1_sct, u1_sct_test, cen_type) {
  
  set.seed(1212)
  if(cen_type ==  "all"){
    data <- u1
    data_test <- u1_test
  } else if(cen_type == "censct"){
    data <- u1_sct
    data_test <- u1_sct_test
  }
  
  stack_surv_test_365 <- brier_365 <- list()
  data_test <- data_test %>% filter(is.na(time) == FALSE)
  
  
  for(xx in 1 : length(tux)){
    t_star <- tux[xx]
    stack_surv_test_365[[xx]] <- rep(NA, nrow(data_test))
    
    for(vv in seq_len(nrow(data_test))){
      grps <- data_test$MOD_RISK[vv]
      grps
      if(grps == "Adverse"){
        rd <- u1_surv$Adverse  
        if(length(which(rd$time <= t_star)) > 0){
          a <- rd$surv[tail(which(rd$time <= t_star), 1)]
        } else {
          a <- rd$surv[1]
        }  
      } else if(grps == "Intermediate"){
        rd <- u1_surv$Intermediate  
        if(length(which(rd$time <= t_star)) > 0){
          a <- rd$surv[tail(which(rd$time <= t_star), 1)]
        } else {
          a <- rd$surv[1]
        }  
      } else if(grps == "Favorable"){
        rd <- u1_surv$Favorable  
        if(length(which(rd$time <= t_star)) > 0){
          a <- rd$surv[tail(which(rd$time <= t_star), 1)]
        } else {
          a <- rd$surv[1]
        }  
      }
      stack_surv_test_365[[xx]][vv] <- a
    }
    
    
    iin <- which(is.na(stack_surv_test_365[[xx]]) == FALSE) 
    
    brier_365[[xx]] <- SurvMetrics::Brier(Surv(data_test$time[iin], data_test$censor[iin]), pre_sp = stack_surv_test_365[[xx]][iin], t_star = t_star)
    brier_365[[xx]]
  }
  # test data res ends 
  
  #ibrier + cumu AUC
  stack_surv_test <- matrix(NA, nrow(data_test), length(tu)) 
  
  for(vv in 1 : nrow(data_test)){
    
    grps <- data_test$MOD_RISK[vv]
    grps
    if(grps == "Adverse"){
      rd <- u1_surv$Adverse  
      
    } else if(grps == "Intermediate"){
      rd <- u1_surv$Intermediate  
      
    } else if(grps == "Favorable"){
      rd <- u1_surv$Favorable  
    } 
    rd <- rd[which(rd$time <= t_max), ]
    stack_surv_test[vv, which(tu %in% rd$time == TRUE)] <- rd$surv # first one is test data for this subject
    stack_surv_test[vv, ]
    for(ll in 1 : 400){ # run loop to get the values carried forward/backward to fill gaps 
      if(sum(is.na(stack_surv_test[vv, ])) > 0){
        if(1 %in% which(is.na(stack_surv_test[vv, ])) ){
          stack_surv_test[vv, 1] <- 1
        }  
        if(1 %!in% which(is.na(stack_surv_test[vv, ])) ){
          stack_surv_test[vv, which(is.na(stack_surv_test[vv, ]))] <- stack_surv_test[vv, (which(is.na(stack_surv_test[vv, ])) - 1)] # Last observation carried forward
        } 
      }
    }
  } # ends for all test   
  sum(is.na(stack_surv_test))
  
  ## replace 0 with a very small number 
  which0 <- function(h) { 
    a1 <- which(h == 0)
    if(length(a1) > 0){
      h[a1] <- 1e-7
    }
    h
  }
  stack_surv_test <- try(do.call(rbind, lapply(seq_len(nrow(stack_surv_test)), function(vv)  which0(h = stack_surv_test[vv, ]))))
  ibrier <- try(SurvMetrics::IBS(Surv(data_test$time, data_test$censor), sp_matrix = stack_surv_test, tu)) 
  
  ## Cumulative case dynamic ROC
  t_survivalROC <- unlist(lapply(seq_len(length(tu)), function(tt) {
    vv <- try(survivalROC(Stime = data_test$time, 
                          status = data_test$censor,
                          marker = - log(stack_surv_test[, which(tu %in% tu[tt])]),
                          predict.time = tu[tt],
                          method = "KM", # nearest neighbour estimation (NNE) & alternative is KM
                          span = 0.25 * nrow(data_test)^(-0.20))$AUC) # span don't need for KM
    
    if(class(vv) %!in% "try-error"){
      vv
    } else{
      NA
    }
  }))
  
  t_risksetROC <- unlist(lapply(seq_len(length(tu)), function(tt) {
    vv <- try(risksetROC::risksetROC(Stime = data_test$time, 
                                     status = data_test$censor,
                                     marker = - log(stack_surv_test[, which(tu %in% tu[tt])]),
                                     predict.time = tu[tt], plot = F,
                                     method = "Cox", # schoenfeld alternatives are Cox-PH
                                     span = 0.25 * nrow(data_test)^(-0.20))$AUC) # span don't need for schoenfeld
    
    if(class(vv) %!in% "try-error"){
      vv
    } else{
      NA
    }
  }))
  
  output = list(t_risksetROC = t_risksetROC, t_survivalROC = t_survivalROC, ibrier  = ibrier, brier_365 = brier_365)
  
}  


# simple survival prediction function based on KM for ELN
simple_survival_prediction_eln_km <- function(tu, tux, t_max, u1_surv, u1, u1_test, u1_sct, u1_sct_test, cen_type) {
  
  set.seed(1212)
  if(cen_type ==  "all"){
    data <- u1
    data_test <- u1_test
  } else if(cen_type == "censct"){
    data <- u1_sct
    data_test <- u1_sct_test
  }
  
  stack_surv_test_365 <- brier_365 <- list()
  data_test <- data_test %>% filter(is.na(time) == FALSE)
  
  for(xx in 1 : length(tux)){
    t_star <- tux[xx]
    stack_surv_test_365[[xx]] <- rep(NA, nrow(data_test))
    
    for(vv in seq_len(nrow(data_test))){
      grps <- data_test$ELN_RISK_GROUP[vv]
      grps
      if(grps == "Adverse"){
        rd <- u1_surv$Adverse  
        #rd$surv[tail(which(rd$time <= t_star), 1)]
        if(length(which(rd$time <= t_star)) > 0){
          a <- rd$surv[tail(which(rd$time <= t_star), 1)]
        } else {
          a <- rd$surv[1]
        }  
      } else if(grps == "Intermediate"){
        rd <- u1_surv$Intermediate  
        if(length(which(rd$time <= t_star)) > 0){
          a <- rd$surv[tail(which(rd$time <= t_star), 1)]
        } else {
          a <- rd$surv[1]
        }  
      } else if(grps == "Favorable"){
        rd <- u1_surv$Favorable  
        if(length(which(rd$time <= t_star)) > 0){
          a <- rd$surv[tail(which(rd$time <= t_star), 1)]
        } else {
          a <- rd$surv[1]
        }  
      }
      stack_surv_test_365[[xx]][vv] <- a
    }
    
    iin <- which(is.na(stack_surv_test_365[[xx]]) == FALSE) # say a group has not been selected in risk stratification model | say it is just adverse + intermediate
    
    brier_365[[xx]] <- SurvMetrics::Brier(Surv(data_test$time[iin], data_test$censor[iin]), 
                                          pre_sp = stack_surv_test_365[[xx]][iin], t_star = t_star)
    brier_365[[xx]]
  }
  # test data res ends 
  
  #ibrier + cumu AUC
  stack_surv_test <- matrix(NA, nrow(data_test), length(tu)) 
  
  for(vv in 1 : nrow(data_test)){
    
    grps <- data_test$ELN_RISK_GROUP[vv]
    grps
    if(grps == "Adverse"){
      rd <- u1_surv$Adverse  
      
    } else if(grps == "Intermediate"){
      rd <- u1_surv$Intermediate  
      
    } else if(grps == "Favorable"){
      rd <- u1_surv$Favorable  
    } 
    rd <- rd[which(rd$time <= t_max), ]
    stack_surv_test[vv, which(tu %in% rd$time == TRUE)] <- rd$surv # first one is test data for this subject
    stack_surv_test[vv, ]
    for(ll in 1 : 400){ # run loop to get the values carried forwrd/backward to fill gaps 
      if(sum(is.na(stack_surv_test[vv, ])) > 0){
        if(1 %in% which(is.na(stack_surv_test[vv, ])) ){
          stack_surv_test[vv, 1] <- 1
        }  
        if(1 %!in% which(is.na(stack_surv_test[vv, ])) ){
          stack_surv_test[vv, which(is.na(stack_surv_test[vv, ]))] <- stack_surv_test[vv, (which(is.na(stack_surv_test[vv, ])) - 1)] # Last observation carried forward
        } 
      }
    }
  } # ends for all test   
  sum(is.na(stack_surv_test))
  
  ## replace 0 with a very small number 
  which0 <- function(h) { 
    a1 <- which(h == 0)
    if(length(a1) > 0){
      h[a1] <- 1e-7
    }
    h
  }
  stack_surv_test <- try(do.call(rbind, lapply(seq_len(nrow(stack_surv_test)), function(vv)  which0(h = stack_surv_test[vv, ]))))
  
  ibrier <- try(SurvMetrics::IBS(Surv(data_test$time, data_test$censor), sp_matrix = stack_surv_test, tu)) 
  
  ## Cumulative case dynamic ROC
  t_survivalROC <- unlist(lapply(seq_len(length(tu)), function(tt) {
    vv <- try(survivalROC(Stime = data_test$time, 
                          status = data_test$censor,
                          marker = - log(stack_surv_test[, which(tu %in% tu[tt])]),
                          predict.time = tu[tt],
                          method = "KM", # nearest neighbour estimation (NNE) & alternative is KM
                          span = 0.25 * nrow(data_test)^(-0.20))$AUC) # span don't need for KM
    
    if(class(vv) %!in% "try-error"){
      vv
    } else{
      NA
    }
    
  }))
  median(t_survivalROC, na.rm = T); quantile(t_survivalROC, 0.05, na.rm = T); quantile(t_survivalROC, 0.95, na.rm = T)
  
  
  
  t_risksetROC <- unlist(lapply(seq_len(length(tu)), function(tt) {
    vv <- try(risksetROC::risksetROC(Stime = data_test$time, 
                                     status = data_test$censor,
                                     marker = - log(stack_surv_test[, which(tu %in% tu[tt])]),
                                     predict.time = tu[tt], plot = F,
                                     method = "Cox", # schoenfeld alternatives are Cox-PH
                                     span = 0.25 * nrow(data_test)^(-0.20))$AUC) # span don't need for schoenfeld
    
    if(class(vv) %!in% "try-error"){
      vv
    } else{
      NA
    }
    
  }))
  median(t_risksetROC, na.rm = T); quantile(t_risksetROC, 0.05, na.rm = T); quantile(t_risksetROC, 0.95, na.rm = T)
  
  
  output = list(t_risksetROC = t_risksetROC, t_survivalROC = t_survivalROC, ibrier  = ibrier, brier_365 = brier_365)
  
}  

