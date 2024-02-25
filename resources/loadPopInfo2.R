masterFile <- "https://docs.google.com/spreadsheets/d/1j5cVIpu5RYRvE1fRHGOv9SliO1JW8-6Dk6l3sdGAHD4/edit?usp=sharing"
master <- gsheet::gsheet2tbl(masterFile)
df <- as.data.frame(master)

if(version$major == "4"){

    df2 <- df[is.na(df$Exclude_reason),]
    
    info <- df2[,c("Common name", "Country","my_locality2", "superpop", "Serial number","depth")]

} else {
    df2 <- df[df$Exclude_reason == "",]

    info <- df2[,c("Common.name", "Country","my_locality2", "superpop", "Serial.number","depth")]

}
colnames(info) <- c("Common_name", "Country", "Locality", "Superpop", "Serial_number", "Depth")

popord <- c("Etosha", "Ovita", "Chobe", "Shangani", "Mana Pools", "Kafue", "Luangwa", "Selous", "Ugalla", "Masai Mara", "Tsavo", "Samburu", "Lake Mburu")
popord <- c("Etosha", "Ovita", "Chobe", "Shangani", "Mana Pools", "Kafue", "Luangwa", "Selous", "Ugalla", "Masai Mara", "Tsavo", "Samburu", "Lake Mburo")



impala_colors <- c("Etosha" = "#8A0E83", "Ovita" = "#F72B69", "Chobe" = "#B81E22", "Shangani" = "#FF927C", "Mana Pools" = "#DB9139", "Kafue" = "#E0D152", "Luangwa" = "#C18E8E", "Selous" = "#698065", "Ugalla" = "#86ABCB", "Masai Mara" = "#26C6C9", "Tsavo" = "#94D9DC", "Samburu" = "#2B597D", "Lake Mburu" = "#0559E8")
impala_colors <- c("Etosha" = "#8A0E83", "Ovita" = "#F72B69", "Chobe" = "#B81E22", "Shangani" = "#FF927C", "Mana Pools" = "#DB9139", "Kafue" = "#E0D152", "Luangwa" = "#C18E8E", "Selous" = "#698065", "Ugalla" = "#86ABCB", "Masai Mara" = "#26C6C9", "Tsavo" = "#94D9DC", "Samburu" = "#2B597D", "Lake Mburo" = "#0559E8")


impala_colors <- impala_colors[popord]




cat("  loaded variables:\n    df: master file\n    info: df with pop info after excluding bad samples \n    impala_colors: vector with colors per locality\n    popord: vector with locality order\n")

