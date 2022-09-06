# 封装的功能性函数============
library('tidyverse')
getname <- function(f_c){ 
  #返回一个向量，内容是所有位置的元素名称，如：
  #  A1   A2   A3   A4    B   X1   X2   X3   X4   X5   X6 
  # "Cs" "Cs" "Cs" "Cs" "Ge" "Br" "Br" "Br" "Br"  "I"  "I" 
  #可以用getname('CsCsCsCsGeBrBrBrBrII')来测试
  sp1 <- strsplit(f_c, split = '')[[1]]
  nameVector <- rep('0', 11)
  id1 = 1
  for (s in 1:length(sp1)) {
    if(s == length(sp1)){
      nameVector[id1] <- sp1[s]
      id1 = 0
    }else if(sp1[s] %in% LETTERS) {
      if(sp1[s + 1] %in% letters){
        nameVector[id1] <- str_c(sp1[s], sp1[s + 1], collapse = '')
        id1 = id1 + 1
      }else{
        nameVector[id1] <- sp1[s]
        id1 = id1 + 1
      }
    }
  }
  names(nameVector) <- c('A1', 'A2', 'A3', 'A4', 'B', 
                         'X1', 'X2', 'X3', 'X4', 'X5', 'X6')
  return(nameVector)
}

getinfo <- function(path){ 
  #输入路径，提取homo，lumo，先判断后提取，如果精度不够则提取为NA
  if("t.outmol" %in% dir(path)){
    outmol <- read_lines(str_c(path, '/t.outmol', collapse = ''))  
  }else{
    return(NA)
  }
  
  if(!("Message: DMol3 job finished successfully" %in% outmol)){
    return(NA) 
  }
  
  lind <- max(1, which("Message: DMol3 job finished successfully" == outmol) - 2000)
  out1 <- outmol[lind: which("Message: DMol3 job finished successfully" == outmol)]
  out2 <- out1
  out3 <- out1
  
  if(length(out2[str_detect(out2, " Yes ")]) == 0){
    return(NA)
  }
  
  #  detect yes or no
  while(length(out3[str_detect(out3, " Yes | No ")]) > 3){
    out3 <- out3[-1]
  }
  
  # extract
  while(length(out2[str_detect(out2, "LUMO")]) > 1){
    out2 <- out2[-1]
  }
  
  if(length(out3[str_detect(out3, " Yes ")]) == 3){
    homo <- as.numeric(gsub('*eV',  '', tail(strsplit(out2[which(str_detect(out2, "LUMO")) - 4], split = ' ')[[1]], 1)))
    lumo <- as.numeric(gsub('*eV',  '', tail(strsplit(out2[which(str_detect(out2, "LUMO")) - 3], split = ' ')[[1]], 1)))
    return(list(homo, lumo))
  }else{
    res1 <- sapply(1:length(out3[str_detect(out3, " No ")]), function(q){
      t1 <- as.double(strsplit(out3[str_detect(out3, " No ")][q], split = "\\|")[[1]])
      t2 <- strsplit(as.character(formatC(t1[!is.na(t1)], format = "e")), split = "e")
      if(as.numeric(t2[[1]][2]) <= as.numeric(t2[[2]][2])){
        return(TRUE)
      }else{
        return(FALSE)
      }
    })
    
    if(sum(res1) == length(res1)){
      homo <- as.numeric(gsub('*eV',  '', tail(strsplit(out2[which(str_detect(out2, "LUMO")) - 4], split = ' ')[[1]], 1)))
      lumo <- as.numeric(gsub('*eV',  '', tail(strsplit(out2[which(str_detect(out2, "LUMO")) - 3], split = ' ')[[1]], 1)))
      return(list(homo, lumo))
    }else{
      return(NA)
    }
  }
}

get_solvent_information <- function(path){
  # 输入路径，获取溶剂效应
  inpuotFile <- read_lines(dir(path, full.names = T)[str_detect(dir(path), '.input')])
  sp1 <- str_split(inpuotFile[str_detect(inpuotFile, 'COSMO_Dielectric')], ' ')[[1]]
  as.numeric(sp1[sp1 != ''][2])
}

# 工作流程============
st1 <- read_csv(file = 'rawdata/A4BX6train_and_test.csv') #获取输入数据

st1$solvent = NA # 在输入数据后面加三列
st1$homo = NA
st1$lumo = NA

d1 <- lapply(1:10, function(q){ #获取所有需要提取的文件的路径
  return(dir(str_c('rawdata/update-0D/', q), full.names = T)  )
})
d1 <- do.call(c ,d1)
ds <- c(dir('rawdata/finish-update-0D', full.names = T), d1)

for(pat in ds[1:length(ds)]){ #遍历路径开始提取
  path <- pat
  name_split <- getname(tail(strsplit(path, '/')[[1]], 1))
  total_energy <- getinfo(path)
  if(is.na(total_energy)){
    total_energy <- list(NA, NA)
  }else{
    print(total_energy)
  }
  
  # 匹配到对应的格子，填入数据
  st1[st1$A1_Symbol == name_split["A1"] & st1$A2_Symbol == name_split["A2"] &
        st1$A3_Symbol == name_split["A3"] & st1$A4_Symbol == name_split["A4"] &
        st1$B_Symbol == name_split["B"] & st1$X1_Symbol == name_split["X1"] &
        st1$X2_Symbol == name_split["X2"] & st1$X3_Symbol == name_split["X3"] &
        st1$X4_Symbol == name_split["X4"] & st1$X5_Symbol == name_split["X5"] &
        st1$X6_Symbol == name_split["X6"], "solvent"] <- get_solvent_information(pat)
  st1[st1$A1_Symbol == name_split["A1"] & st1$A2_Symbol == name_split["A2"] &
        st1$A3_Symbol == name_split["A3"] & st1$A4_Symbol == name_split["A4"] &
        st1$B_Symbol == name_split["B"] & st1$X1_Symbol == name_split["X1"] &
        st1$X2_Symbol == name_split["X2"] & st1$X3_Symbol == name_split["X3"] &
        st1$X4_Symbol == name_split["X4"] & st1$X5_Symbol == name_split["X5"] &
        st1$X6_Symbol == name_split["X6"], "homo"] <- total_energy[[1]][1]
  st1[st1$A1_Symbol == name_split["A1"] & st1$A2_Symbol == name_split["A2"] &
        st1$A3_Symbol == name_split["A3"] & st1$A4_Symbol == name_split["A4"] &
        st1$B_Symbol == name_split["B"] & st1$X1_Symbol == name_split["X1"] &
        st1$X2_Symbol == name_split["X2"] & st1$X3_Symbol == name_split["X3"] &
        st1$X4_Symbol == name_split["X4"] & st1$X5_Symbol == name_split["X5"] &
        st1$X6_Symbol == name_split["X6"], "lumo"] <- total_energy[[2]][1]
}

write_csv(st1, file = 'rawdata/A4BX6train_and_test_fill1.csv')
