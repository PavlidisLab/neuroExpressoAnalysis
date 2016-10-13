library(XLConnect)


download.file("http://www.cell.com/cms/attachment/2021775604/2041641016/mmc1.zip",
              destfile='data-raw/Mouse_Cell_Type_Data//s100a10AllowedProbes.zip')
system('unzip data-raw/Mouse_Cell_Type_Data/s100a10AllowedProbes.zip -d data-raw/Mouse_Cell_Type_Data')


probes = loadWorkbook('data-raw/Mouse_Cell_Type_Data/cell_6244_mmc1.xlsx')
probes = readWorksheet(probes, sheet = 1, header = TRUE)

allowedProbesS100a10 = probes$Probe.Set.ID %>% gemmaGeneMatch(chipFile='data-raw/GemmaAnnots/GPL1261') %>% unique

devtools::use_data(allowedProbesS100a10, overwrite = TRUE)