tax <- tibble::tribble(
  ~strainID,             ~genus,        ~species,
  "HAMBI_0006",      "Pseudomonas",        "putida",
  "HAMBI_0097",    "Acinetobacter",       "lwoffii",
  "HAMBI_0105",    "Agrobacterium",   "tumefaciens",
  "HAMBI_0262",    "Brevundimonas",       "bullata",
  "HAMBI_0403",        "Comamonas",  "testosteroni",
  "HAMBI_1279",           "Hafnia",         "alvei",
  "HAMBI_1287",      "Citrobacter",        "koseri",
  "HAMBI_1292",       "Morganella",      "morganii",
  "HAMBI_1299",         "Kluyvera",    "intermedia",
  "HAMBI_1842",      "Sphingobium",    "yanoikuyae",
  "HAMBI_1896", "Sphingobacterium",  "spiritivorum",
  "HAMBI_1923", "Myroides",              "odoratus",
  "HAMBI_1972",        "Aeromonas",        "caviae",
  "HAMBI_1977",      "Pseudomonas",  "chlororaphis",
  "HAMBI_1988",     "Chitinophaga",        "sancti",
  "HAMBI_2159", "Paraburkholderia",   "caryophylli",
  "HAMBI_2160",       "Bordetella",         "avium",
  "HAMBI_2164",      "Cupriavidus",       "necator",
  "HAMBI_2443",       "Paracoccus", "denitrificans",
  "HAMBI_2494", "Paraburkholderia",   "kururiensis",
  "HAMBI_2659", "Stenotrophomonas",   "maltophilia",
  "HAMBI_2792",        "Moraxella",         "canis",
  "HAMBI_3031",         "Niabella",  "yanshanensis",
  "HAMBI_3237",       "Microvirga",   "lotononidis"
)

rRNA_counts <- tibble::tribble(
                ~strainID, ~copies,
             "HAMBI_0006",               7L,
             "HAMBI_0097",               7L,
             "HAMBI_0105",               4L,
             "HAMBI_0262",               3L,
             "HAMBI_0403",               9L,
             "HAMBI_1279",               7L,
             "HAMBI_1287",               7L,
             "HAMBI_1292",               7L,
             "HAMBI_1299",               8L,
             "HAMBI_1842",               4L,
             "HAMBI_1896",               4L,
             "HAMBI_1923",               6L,
             "HAMBI_1972",               10L,
             "HAMBI_1977",               5L,
             "HAMBI_1988",               5L,
             "HAMBI_2159",               4L,
             "HAMBI_2160",               3L,
             "HAMBI_2164",               5L,
             "HAMBI_2443",               3L,
             "HAMBI_2494",               4L,
             "HAMBI_2659",               4L,
             "HAMBI_2792",               4L,
             "HAMBI_3031",               2L,
             "HAMBI_3237",               6L
             )

hambi_colors <- c(`HAMBI_1287` = "#6DA3EE", `HAMBI_1977` = "#F1C82E", `HAMBI_0403` = "#FB8A5C", 
                  `HAMBI_2659` = "#624090", `HAMBI_1972` = "#26818E", `HAMBI_1292` = "#32702C", 
                  `HAMBI_1923` = "#FC2EDB", `HAMBI_1279` = "#FE1C35", `HAMBI_1299` = "#BEEF60", 
                  `HAMBI_1896` = "#1C26FB", `HAMBI_0097` = "#F0DACB", `HAMBI_1988` = "#D71C76", 
                  `HAMBI_0006` = "#870DAE", `HAMBI_2792` = "#0DFE32", `HAMBI_3031` = "#FC94D1", 
                  `HAMBI_2160` = "#16EEA5", `HAMBI_0105` = "#00F5F7", `HAMBI_3237` = "#956616", 
                  `HAMBI_2494` = "#DABBF3", `HAMBI_2164` = "#DF26FD", `HAMBI_0262` = "#8D324F", 
                  `HAMBI_2443` = "#9DE0C5", `HAMBI_2159` = "#92950D", `HAMBI_1842` = "#F489FD"
)

# opposite of %in% fuction
`%nin%` = Negate(`%in%`)

# logit transform
logit = function(x){
  log(x/(1-x))
}

minnz = function(V) {
    # Calculates the smallest value of the vector except for 0 (non-zero minumum)
    # Argument: vector
    C <- NULL        # prepare
    k <- length(V)   # count to
    for (i in 1:k) { # check all
      if ((V[i] == 0) == FALSE) (C[i] <- V[i]) else (C[i] <- 9999919) # if V[i] is not 0, add it to C
    }
    m <- min(C)               # minimum of V, not counting 0
    if (max(V) == 1) (m <- 1) # fix for binary vectors (0,1)
    if (m == 9999919) (warning("Error: Minimum calculation failed."))  # warning because of hard-coded replacement
    return(m)
  }

quibble95 = function(x, q = c(0.025, 0.5, 0.975)) {
  tibble(x = quantile(x, q), quantile = c("q2.5", "q50", "q97.5"))
}

`%nin%` <- Negate(`%in%`)

getfiles <- function(dir, pattern){
  list.files(dir, full.names = TRUE, pattern = pattern)
}

getnames <- function(fileslist, myregex){
  set_names(fileslist,
            str_extract(
              fileslist,
              regex(myregex)
            ))
}

# "(?<=[/])([^/]+)(?=(_[:alpha:]+_[:alpha:]+\\.)[^.]+)|(?<=[/])([^/]+)(?=(_[:alpha:]+\\.)[^.]+)"
