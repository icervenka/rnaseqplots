mito_atp_replacer <- c(
  "Atp5f1a" = "Atp5a1",
  "Atp5f1b" = "Atp5b",
  "Atp5f1c" = "Atp5c1",
  "Atp5f1d" = "Atp5d",
  "Atp5mc1" = "Atp5g1",
  "Atp5mc2" = "Atp5g2",
  "Atp5mc3" = "Atp5g3",
  "Atp5mf" = "Atp5j2",
  "Atp5mg" = "Atp5l",
  "Atp5pd" = "Atp5h",
  "Atp5pf" = "Atp5j",
  "Atp5po" = "Atp5o",
  "Cyb" = "Cytb"
)

usethis::use_data(mito_atp_replacer, overwrite = TRUE)