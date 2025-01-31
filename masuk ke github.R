# Set working directory
setwd("D:/Kuliah S2/Thesis/Syntax")

#Load R package
library(GWmodel)      ## GW models
library(plyr)         ## Data management
library(sp)           ## Spatial Data management
library(spdep)        ## Spatial autocorrelation
library(RColorBrewer) ## Visualization
library(classInt)     ## Class intervals
library(raster)       ## spatial data
library(grid)         ## plot
library(gridExtra)    ## Multiple plot
library(ggplot2)      #  plotting
library(tidyverse)    # data 
library(SpatialML)    # Geographically weigted regression
library(dplyr)
library(sf)

# Define data folder
dataFolder<-"D:/Kuliah S2/Thesis/Syntax"
county<-shapefile(paste0("D:/Kuliah S2/Thesis/Syntax/[LapakGIS.com] Batas Wilayah Kabupaten 2024/LapakGIS_Batas_Kabupaten_2024.shp"))
jawa_provinces <- c("Banten", "DKI Jakarta", "Jawa Barat", "Jawa Tengah", "Daerah Istimewa Yogyakarta", "Jawa Timur")
state <- list("sp.lines", as(county[county$WADMPR %in% jawa_provinces, ], "SpatialLines"), col="grey50", lwd=.5,lty=3) 
df<-read.csv(paste0("D:/Kuliah S2/Thesis/Syntax/Data AKB1.csv"), header=T, sep=";")
str(df)


#===== Penggabungan Data dan SHP
county1 <- st_read("D:/Kuliah S2/Thesis/Syntax/[LapakGIS.com] Batas Wilayah Kabupaten 2024/LapakGIS_Batas_Kabupaten_2024.shp")
county1 <- st_make_valid(county1)  # Validasi data spasial
kab.kot_jawa <- df$WADMKK
county_jawa <- county1 %>% filter(WADMKK %in% df$WADMKK)
county_jawa <- county1 %>% filter(WADMKK %in% kab.kot_jawa)

# Gabungkan data AKB dengan shapefile berdasarkan kolom WADMKK
data_sf <- county_jawa %>%
  left_join(df, by = "WADMKK") %>%
  select(-any_of(c("METADATA", "UPDATED")))
head(data_sf)

#=== Filter Tahun 2022
data22 <- data_sf %>%  
  filter(Tahun == 2022)
data22 <- data22[,c(3:19,29)] #variabel yang digunakan
head(data22)

#========= EKSPLORASI DATA
#==== DESKRIPSI DATA
library(psych)
des <- describe(st_drop_geometry(data22))
#==== KORELASI
korelasi <- cor(st_drop_geometry(data22[,6:17]))
korelasi

library(ggcorrplot)
png("22_Korelasi1.png")
ggcorrplot(korelasi, lab = TRUE, colors = c("red", "white", "blue"))
dev.off()

#=== AUTOKORELASI SPASIAL
library(spdep)
coords22 <- cbind(data22$x, data22$y) 
nb <- knn2nb(knearneigh(coords22, k = 5))  
listw <- nb2listw(nb, style = "W")
morans_test <- moran.test(data22$AKB, listw)
print(morans_test)

#============== SCALE DATA VARIABEL PREDIKTOR
data_22 <- st_drop_geometry(data22) 
data_22[, c(7:17)] = scale(data_22[, c(7:17)]);head(data_22) 
coords_sf22 <- st_as_sf(data_22, coords = c("x", "y"), crs = 4326) #data point
print(coords_sf22)

write_xlsx(data_22, path = "Scale Data Prediktor.xlsx")

#============== SPATIAL STRATIFIDIED SAMPLING

#============== CROSS-VALIDATION 2022
#==== Blok Cross Validation
library(blockCV)
set.seed(42)
cv_folds22 <- cv_spatial(
  x = coords_sf22,          # Data spasial
  k = 5,                    # Jumlah fold yang diinginkan
  selection = "systematic", # Pemilihan sistematis fold
  hexagon = T,              # Blok hexagonal
  rows_cols = c(53, 53),    # Baris dan kolom untuk pembagian blok
  iteration = 100,          # Banyaknya iterasi untuk pemilihan fold acak
  extend = 0.5,             # Ekstensi blok
  progress = TRUE,          # Menampilkan progress bar
  report = TRUE,            # Menampilkan laporan fold
  plot = TRUE,              # Menampilkan plot fold
  offset = c(0.1, 0.1)      # Geser sedikit untuk mendistribusikan fold lebih baik
)
cv_folds22$records

# Assign fold IDs to the data_22 dataframe
data22$fold_id <- cv_folds22$folds_ids #data point
data_22$fold_id <- cv_folds22$folds_ids #data multipoligon
head(data22)

unique_folds <- unique(data22$fold_id)

#=== SPATIAL CV
# Initialize list to store results
fold_results <- list()

for (fold in unique_folds) {
  print(paste("Processing fold:", fold))
  
  # Data train and test
  train_data <- data_22[data_22$fold_id != fold, ]
  test_data <- data_22[data_22$fold_id == fold, ]
  
  print(paste("Train data size for fold", fold, ":", nrow(train_data)))
  print(paste("Test data size for fold", fold, ":", nrow(test_data)))
  
  if (nrow(train_data) == 0 || nrow(test_data) == 0) {
    print(paste("Fold", fold, "has insufficient data for training/testing"))
    next  # Skip the fold if it has no sufficient data
  }
  
  # Define target and predictors
  target <- "AKB"
  predictors <- c("PK", "RLS", "PKKP", "KP", "SL", "AML", 
                  "IDL", "ASI", "BBLR", "Rdokter", "Rpuskesmas")
  
  # Convert to data frame with relevant columns
  train_data_df <- data.frame(train_data[, c(target, predictors)])  # Training data
  test_data_df  <- data.frame(test_data[, c(target, predictors)])   # Testing data
  
  # Extract coordinates from the original data
  train_coords <- as.matrix(train_data[, c("x", "y")])  # Training coordinates
  test_coords  <- as.matrix(test_data[, c("x", "y")])   # Testing coordinates
  
  # Add coordinates as 'x' and 'y' columns to the training and testing data frames
  train_data_df$x <- train_coords[, 1]
  train_data_df$y <- train_coords[, 2]
  test_data_df$x <- test_coords[, 1]
  test_data_df$y <- test_coords[, 2]
  
  # Save data frames in the list for further inspection
  fold_results[[paste0("fold_", fold)]] <- list(train = train_data_df, test = test_data_df)
  
  # Now that the data is properly set, proceed with bandwidth optimization
  print(paste("Optimizing bandwidth for fold", fold))
  bw_result <- tryCatch({
    grf.bw(
      formula = paste("AKB ~", paste(predictors, collapse = " + ")),
      dataset = train_data,           # Dataset for the model
      kernel = "adaptive",               # Kernel type
      coords = train_data_df[, c("x", "y")],  # Spatial coordinates
      bw.min = 15,                       # Minimum bandwidth
      bw.max = nrow(train_data)/3,    # Maximum bandwidth
      step = 1,                          # Step size for bandwidth iteration
      trees = 50,                        # Number of trees in the local forest
      mtry = 4,                          # Number of variables for split
      importance = "permutation",        # Mode for variable importance
      nthreads = parallel::detectCores(),  # Number of threads for parallel computation
      forests = FALSE,                   # Do not store all local forests
      geo.weighted = TRUE                # Use spatial weighting
    )
  }, error = function(e) {
    print(paste("Error in bandwidth optimization for fold", fold, ":", e$message))
    return(NULL)
  })
  
  if (is.null(bw_result)) {
    print(paste("Skipping fold", fold, "due to error in bandwidth optimization"))
    next  # Skip the fold if bandwidth optimization fails
  }
  
  # Use the optimal bandwidth value from bw_result
  optimal_bw <- bw_result$Best.BW
  print(paste("Optimal bandwidth for fold", fold, ":", optimal_bw))
  
  # 7. Perform Random Search for optimal ntree
  print(paste("Performing Random Search for ntree for fold", fold))
  optimal_ntree <- random_search_ntree(train_data_df, target, predictors, optimal_bw, train_data_df[, c("x", "y")])
  print(paste("Optimal ntree for fold", fold, ":", optimal_ntree))
  
  # 8. Build model using grf with optimized bandwidth and ntree
  print(paste("Building model for fold", fold))
  model <- tryCatch({
    grf(
      formula = as.formula(paste(target, "~", paste(predictors, collapse = "+"))),
      dframe = train_data_df,   # Training data in data.frame format
      bw = optimal_bw,          # Bandwidth from the optimization result
      kernel = "adaptive",      # Kernel type
      coords = train_data_df[, c("x", "y")],  # Coordinates
      mtry = 4,
      importance = "permutation",  
      ntree = optimal_ntree,    # Optimized ntree from Random Search
      nthreads = parallel::detectCores(), # Jumlah thread untuk paralelisme
      forests = TRUE,                    # Tidak menyimpan semua local forest
      geo.weighted = T,                # Menggunakan pembobotan spasial
      print.results = TRUE,                # Menampilkan hasil analisis
      seed = 1234567
    )
  }, error = function(e) {
    print(paste("Error in model training for fold", fold, ":", e$message))
    return(NULL)
  })
  
  if (is.null(model)) {
    print(paste("Skipping fold", fold, "due to error in model training"))
    next  # Skip the fold if model training fails
  }
  
  print(paste("Model successfully trained for fold:", fold))
  
  # 9. Make predictions on test data
  print(paste("Making predictions for fold", fold))
  predictions <- tryCatch({
    predict.grf(model, 
                new.data = test_data_df[, c(target, predictors, "x", "y")],
                x.var.name = "x",   # Column name for X coordinates
                y.var.name = "y",   # Column name for Y coordinates
                local.w = 1,        # Weight for local model
                global.w = 0        # Weight for global model
    )
  }, error = function(e) {
    print(paste("Error in prediction for fold", fold, ":", e$message))
    return(NA)
  })
  
  # Skip fold if predictions contain NA
  if (any(is.na(predictions))) {
    print(paste("Skipping fold", fold, "due to NA values in predictions"))
    next
  }
  
  # 10. Evaluate model performance
  rmse <- sqrt(mean((test_data_df[[target]] - predictions)^2))
  rsq <- 1 - sum((test_data_df[[target]] - predictions)^2) / 
    sum((test_data_df[[target]] - mean(train_data_df[[target]]))^2)
  
  # 11. Output evaluation for each fold
  print(paste("RMSE for fold", fold, ":", rmse))
  print(paste("R-squared for fold", fold, ":", rsq))
  
  # 12. Get local predictions and MSE for training data from model
  train_rsq <- model[["LocalModelSummary"]][["l.r.Pred"]]  # Local R-squared for training data
  train_rmse <- sqrt(model[["LocalModelSummary"]][["l.MSE.Pred"]])  # RMSE from training data
  train_roob <- model[["LocalModelSummary"]][["l.r.OOB"]]
  
  # 14. Return the result for the current fold, including training metrics
  print(paste("Training RMSE for fold", fold, ":", train_rmse))
  print(paste("Training R-squared for fold", fold, ":", mean(train_rsq)))
  
  # Return the result for the fold
  fold_results[[paste0("fold_", fold)]]$metrics <- list(
    rmse_test = rmse,
    rsq_test = rsq,
    rmse_train = train_rmse,
    rsq_train = mean(train_rsq),
    roob_train = train_roob,
    opt_ntree= optimal_ntree,
    opt_bw = optimal_bw
  )
}

hasil22 <- rbind(
  data.frame(fold = 1, fold_results$fold_1$metrics),
  data.frame(fold = 2, fold_results$fold_2$metrics),
  data.frame(fold = 3, fold_results$fold_3$metrics),
  data.frame(fold = 4, fold_results$fold_4$metrics),
  data.frame(fold = 5, fold_results$fold_5$metrics)
)
hasil22


############################################### GAK DIMASUKIN KE SCRIPT TESIS
#============= HASIL SPLITING
#==== Mengambil nilai fold dari baris dengan R² test tertinggi
best_fold22 <- hasil22[which.max(hasil22$rsq_test), "fold"]
print(best_fold22)

#==== data point ==> untuk pemodelan
train.df22 <- data_22[data_22$fold_id != best_fold22, ];#head(train.df22)  # Training data
test.df22  <- data_22[data_22$fold_id == best_fold22, ]; #head(test.df22) # Testing data

#==== data multipolygon
train22 <- data22[data22$fold_id != best_fold22, ];head(train22)  # Training data
test22  <- data22[data22$fold_id == best_fold22, ]; head(test22) # Testing data

############################################### GAK DIMASUKIN KE SCRIPT TESIS

#==== Menambahkan kolom Tipe Data Training/Testing
data22$Data_Type <- ifelse(data22$WADMKK %in% train.df22$WADMKK, "Training", 
                           ifelse(data22$WADMKK %in% test.df22$WADMKK, "Testing", NA))
head(data22)
#==== Deskripsi Data Training dan Testing
des2 <- describe(st_drop_geometry(data22[data22$Data_Type == "Training", ]))
des3 <- describe(st_drop_geometry(data22[data22$Data_Type == "Testing", ]))

# Membuat peta
g22 <-ggplot() +
  geom_sf(data = train22, aes(color = "Training"), size = 2, alpha = 0.7) +  # Peta data training
  geom_sf(data = test22, aes(color = "Testing"), size = 2, alpha = 0.7) +    # Peta data testing
  scale_color_manual(values = c("Training" = "blue", "Testing" = "red")) +   # Warna untuk training dan testing
  labs(title = "Peta Distribusi Data Training dan Testing",
       color = "Kategori") +
  theme_minimal()
ggsave("D:/Kuliah S2/Thesis/Syntax/Eksplorasi/22__train test.png", plot =g22, dpi = 600, width = 8, height = 6, units = "in")


#================= GEOGRAPHICALLY WEIGHTED RANDOM FOREST
#====== BANDWIDTH OPTIMUM
set.seed(42)
target <- "AKB"
predictors <- c("PK", "RLS", "PKKP", "KP", "SL", "AML", 
                "IDL", "ASI", "BBLR", "Rdokter", "Rpuskesmas")

bw_result22 <- grf.bw(
  formula = as.formula(paste(target, "~", paste(predictors, collapse = "+"))),
  dataset = train.df22,           # Dataset untuk model
  kernel = "adaptive",               # Jenis kernel
  coords = train.df22[, c("x", "y")],  # Matriks koordinat spasial (X,Y)
  bw.min = 15,                       # Minimum bandwidth
  bw.max = nrow(train.df22)/3,    # Maximum bandwidth
  step = 1,                          # Langkah iterasi bandwidth
  trees = 50,                        # Jumlah pohon untuk setiap local forest
  mtry = 4,                          # Jumlah variabel yang dipilih secara acak untuk split
  importance = "permutation",        # Mode penghitungan pentingnya variabel
  nthreads = parallel::detectCores(),  # Jumlah thread untuk paralelisme
  forests = FALSE,                   # Tidak menyimpan semua local forest
  geo.weighted = TRUE                # Menggunakan pembobotan spasial
)

#==== Hasil Bandwidth
bw.opt22 <- data.frame(bw_result22$tested.bandwidths)

# 10 terbaik
bw.opt22_sorted <- bw.opt22[order(-bw.opt22$Local), ]
top_10_bw <- head(bw.opt22_sorted, 10)
print(top_10_bw)

#==== Bandwidth terbaik
optimal_bw22 <- bw_result22$Best.BW
optimal_bw22

#==== Ploting Bandwidth
bw <- bw.opt22
max_value <- bw$Local[bw$Bandwidth == optimal_bw22]
optimal_bandwidth <- bw$Bandwidth[bw$Bandwidth == optimal_bw22]

band.opt <- ggplot(bw, aes(x = Bandwidth)) +
  geom_line(aes(y = Local, color = "Local"), size = 1.2) +  # Garis lebih tebal
  annotate("point", x = optimal_bandwidth, y = max_value, color = "red", size = 3) +  # Titik maksimum
  annotate(
    "text",
    x = optimal_bandwidth, y = max_value,
    label = paste0("Bandwidth: ", optimal_bandwidth),
    vjust = -1, hjust = 0.5, color = "black", size = 4  # Penyesuaian posisi teks
  ) +
  labs(
    title = "Kinerja Bandwidth",
    x = "Bandwidth",
    y = "R-Square"
  ) +
  theme_minimal() +
  theme(legend.position = "none")  # Menghilangkan legenda

ggsave("D:/Kuliah S2/Thesis/Syntax/Eksplorasi/22_Bandwidth Optimum.png", plot =band.opt, dpi = 600, width = 8, height = 6, units = "in")

#==== MODEL GWRF
optimal_ntree22 <- hasil22[which.max(hasil22$rsq_test), "opt_ntree"]
optimal_ntree22

model.gwrf22 <- grf(
  formula = as.formula(paste(target, "~", paste(predictors, collapse = "+"))),
  dframe = train.df22,   # Training data in data.frame format
  bw = optimal_bw22,          # Bandwidth from the optimization result
  kernel = "adaptive",      # Kernel type
  coords = train.df22[, c("x", "y")],  # Coordinates
  mtry = 4,
  importance = "permutation",  
  ntree = optimal_ntree22,    # Optimized ntree from Random Search
  nthreads = parallel::detectCores(), # Jumlah thread untuk paralelisme
  forests = TRUE,                    # Tidak menyimpan semua local forest
  geo.weighted = T,                # Menggunakan pembobotan spasial
  print.results = TRUE,                # Menampilkan hasil analisis
  seed = 1234567
)
summary(model.gwrf22$Global.Model)

#========= POHON KEPUTUSAN GWRF
tree.info <- treeInfo(model.gwrf22[["Forests"]][[2]] , tree=1)

library(DiagrammeR)

write_dot <- function(tree_info, file_name) {  
  # Start the DOT graph  
  cat("digraph Tree {\n", file = file_name)  
  cat("  node [shape = rectangle, style = filled, fillcolor = lightblue]\n", file = file_name, append = TRUE)  
  
  # Loop through each row to create nodes and edges  
  for (i in seq_len(nrow(tree_info))) {  
    # Check if the node is terminal  
    if (tree_info$terminal[i]) {  
      # Label for terminal nodes  
      label <- sprintf("%.2f", tree_info$prediction[i])  
    } else {  
      # Label for non-terminal nodes  
      label <- sprintf("%s < %.2f", tree_info$splitvarName[i], tree_info$splitval[i])  
    }  
    
    # Generate the node with label  
    cat(sprintf("  %d [label=\"%s\"];\n", tree_info$nodeID[i], label), file = file_name, append = TRUE)  
    
    # Generate edges for children if they exist  
    if (!is.na(tree_info$leftChild[i])) {  
      cat(sprintf("  %d -> %d [label = 'Yes'];\n", tree_info$nodeID[i], tree_info$leftChild[i]), file = file_name, append = TRUE)  
    }  
    if (!is.na(tree_info$rightChild[i])) {  
      cat(sprintf("  %d -> %d [label = 'No'];\n", tree_info$nodeID[i], tree_info$rightChild[i]), file = file_name, append = TRUE)  
    }  
  }  
  
  # Close the DOT graph  
  cat("}\n", file = file_name, append = TRUE)  
}  

# Save to .dot file  
write_dot(tree.info, "tree.dot")  

# Visualize using DiagrammeR or in Graphviz Viewer  
grViz("tree.dot")

#==== PREDIKSI
FIPS22.train.xy <- st_as_sf(train22[c(1:6)])
FIPS22.test.xy <- st_as_sf(test22[c(1:6)])

FIPS22.train.xy$Pred.GWRF.Local <- model.gwrf22[["LGofFit"]][["LM_yfitPred"]] 

Pred.GWRF.Local25 <- predict.grf(model.gwrf22 ,test.df22, x.var.name = "x", y.var.name = "y", local.w = 0.25, global.w = 0.75) 
Pred.GWRF.Local5 <- predict.grf(model.gwrf22 ,test.df22, x.var.name = "x", y.var.name = "y", local.w = 0.5, global.w = 0.5) 
Pred.GWRF.Local75 <- predict.grf(model.gwrf22 ,test.df22, x.var.name = "x", y.var.name = "y", local.w = 0.75, global.w = 0.25) 
Pred.GWRF.Local1 <- predict.grf(model.gwrf22 ,test.df22, x.var.name = "x", y.var.name = "y", local.w = 1, global.w = 0) 

# Memilih Bobot Parameter untuk Prediksi
results <- data.frame(
  Model = c('w=0.25', 'w=0.50', 'w=0.75', 'w=1'),
  RMSE = c(
    round(sqrt(mean((Pred.GWRF.Local25 - FIPS22.test.xy$AKB)^2, na.rm = TRUE)), digits = 3),
    round(sqrt(mean((Pred.GWRF.Local5 - FIPS22.test.xy$AKB)^2, na.rm = TRUE)), digits = 3),
    round(sqrt(mean((Pred.GWRF.Local75 - FIPS22.test.xy$AKB)^2, na.rm = TRUE)), digits = 3),
    round(sqrt(mean((Pred.GWRF.Local1 - FIPS22.test.xy$AKB)^2, na.rm = TRUE)), digits = 3)
  ),
  MAE = c(
    round(mean(abs(Pred.GWRF.Local25 - FIPS22.test.xy$AKB), na.rm = TRUE), digits = 3),
    round(mean(abs(Pred.GWRF.Local5 - FIPS22.test.xy$AKB), na.rm = TRUE), digits = 3),
    round(mean(abs(Pred.GWRF.Local75 - FIPS22.test.xy$AKB), na.rm = TRUE), digits = 3),
    round(mean(abs(Pred.GWRF.Local1 - FIPS22.test.xy$AKB), na.rm = TRUE), digits = 3)
  ),
  R2 = c(
    round(summary(lm(FIPS22.test.xy$AKB ~ Pred.GWRF.Local25))$r.squared, digits = 3),
    round(summary(lm(FIPS22.test.xy$AKB ~ Pred.GWRF.Local5))$r.squared, digits = 3),
    round(summary(lm(FIPS22.test.xy$AKB ~ Pred.GWRF.Local75))$r.squared, digits = 3),
    round(summary(lm(FIPS22.test.xy$AKB ~ Pred.GWRF.Local1))$r.squared, digits = 3)
  )
)

results

FIPS22.test.xy$Pred.GWRF.Local <- predict.grf(model.gwrf22 ,test.df22, x.var.name = "x", y.var.name = "y", local.w = 0.5, global.w = 0.5) 

cat('GWRF Local Test RMSE:', round(sqrt(mean((FIPS22.test.xy$Pred.GWRF.Local-FIPS22.test.xy$AKB)^2 , na.rm = TRUE)), digits=3), '\n')
cat('GWRF Local Test MAE:', round(mean(abs(FIPS22.test.xy$Pred.GWRF.Local-FIPS22.test.xy$AKB) , na.rm = TRUE ), digits=3), '\n')
cat('GWRF Local Test R2:', round(summary(lm(AKB~Pred.GWRF.Local5,FIPS22.test.xy))$r.squared, digits=3), '\n')


FIPS22.train.xy$Res.GWRF <- FIPS22.train.xy$AKB-FIPS22.train.xy$Pred.GWRF.Local
FIPS22.test.xy$Res.GWRF <- FIPS22.test.xy$AKB-FIPS22.test.xy$Pred.GWRF.Local
#==== R2 local testing
y_actual_test <- FIPS22.test.xy$AKB
y_pred_test <- FIPS22.test.xy$Pred.GWRF.Local
y_mean_test <- mean(y_actual_test)
R_squared_local <- 1 - ((y_actual_test - y_pred_test)^2 / (y_actual_test - y_mean_test)^2)
R_squared_local
summary(R_squared_local)

#==== R2 local training
y_actual_train <- FIPS22.train.xy$AKB
y_pred_train <- FIPS22.train.xy$Pred.GWRF.Local
y_mean_train <- mean(y_actual_train)
R_squared_local1 <- 1 - ((y_actual_train - y_pred_train)^2 / (y_actual_train - y_mean_train)^2)
R_squared_local1
summary(R_squared_local1)

#==== Menambahkan R2 GWRF ke FIPS22
FIPS22.test.xy$R2.GWRF <-  R_squared_local
FIPS22.train.xy$R2.GWRF <-  R_squared_local1

# Hitung nilai R2 dan RMSE untuk dimasukan ke dalam plot
local_rmse <- round(sqrt(mean((FIPS22.test.xy$Pred.GWRF.Local - FIPS22.test.xy$AKB)^2, na.rm = TRUE)), digits = 3)
local_r2 <- round(summary(lm(AKB ~ Pred.GWRF.Local5, FIPS22.test.xy))$r.squared, digits = 3)
local_mae <- round(mean(abs(FIPS22.test.xy$Pred.GWRF.Local-FIPS22.test.xy$AKB)), digits = 3)

#=== Ploting actual vs prediction Testing
p.test.GWRF22 <- ggplot(data = FIPS22.test.xy, aes(x = AKB, y = Pred.GWRF.Local5)) + 
  geom_point(size = 2.0, color = "blue") +  # Titik berwarna biru
  geom_smooth(method = "lm", se = FALSE, colour = "black") +  # Garis regresi berwarna hitam
  labs(x = "Observed", y = "Predicted", 
       title = "Observed vs Predicted (GWRF Local)") +
  geom_label(aes(
    x = min(FIPS22.test.xy$AKB, na.rm = TRUE) + 0.1 * diff(range(FIPS22.test.xy$AKB, na.rm = TRUE)), 
    y = max(FIPS22.test.xy$Pred.GWRF.Local, na.rm = TRUE) - 0.1 * diff(range(FIPS22.test.xy$Pred.GWRF.Local, na.rm = TRUE)), 
    label = paste0("MAE: ", local_mae, "\nRMSE: ", local_rmse, "\nR²: ", local_r2)),
    hjust = 0.5, vjust = 0.2, size = 4.5, fill = "white", color = "black", label.size = 0.1) +  # Warna kotak putih
  theme_minimal() +  # Latar belakang putih dengan tema minimal
  theme(legend.position = "top")  # Legend di atas

# Menampilkan plot
print(p.test.GWRF22)

ggsave("D:/Kuliah S2/Thesis/Syntax/Eksplorasi/22_GWRF_Actual vs Prediction.png", plot = p.test.GWRF22, dpi = 800,width = 5.2, height = 5.2,  units = "in")

#=== Ploting actual vs prediction Training
traingwrf_rmse <- round(sqrt(mean((FIPS22.train.xy$Pred.GWRF.Local - FIPS22.train.xy$AKB)^2, na.rm = TRUE)), digits = 4)
traingwrf_r2 <- round(summary(lm(AKB ~ Pred.GWRF.Local, FIPS22.train.xy))$r.squared, digits = 4)
traingwrf_mae <- round(mean(abs(FIPS22.train.xy$Pred.GWRF.Local-FIPS22.train.xy$AKB)), digits = 4)

p.test.GWRF22train <- ggplot(data = FIPS22.train.xy, aes(x = AKB, y = Pred.GWRF.Local)) + 
  geom_point(size = 2.0, color = "blue") +  # Titik berwarna biru
  geom_smooth(method = "lm", se = FALSE, colour = "black") +  # Garis regresi berwarna hitam
  labs(x = "Observed", y = "Predicted", 
       title = "Observed vs Predicted (GWRF Local)") +
  geom_label(aes(
    x = min(FIPS22.train.xy$AKB, na.rm = TRUE) + 0.1 * diff(range(FIPS22.train.xy$AKB, na.rm = TRUE)), 
    y = max(FIPS22.train.xy$Pred.GWRF.Local, na.rm = TRUE) - 0.1 * diff(range(FIPS22.train.xy$Pred.GWRF.Local, na.rm = TRUE)), 
    label = paste0("MAE: ", traingwrf_mae, "\nRMSE: ", traingwrf_rmse, "\nR²: ", traingwrf_r2)),
    hjust = 0.5, vjust = 0.2, size = 4.5, fill = "white", color = "black", label.size = 0.1) +  # Warna kotak putih
  theme_minimal() +  # Latar belakang putih dengan tema minimal
  theme(legend.position = "top")  # Legend di atas

# Menampilkan plot
print(p.test.GWRF22train)
ggsave("D:/Kuliah S2/Thesis/Syntax/Eksplorasi/22_GWRF_Actual vs Prediction.png", plot = p.test.GWRF22train, dpi = 800,width = 5.2, height = 5.2,  units = "in")

#==== VIAUALISASI HASIL PREDIKSI
# Gabungkan data training dan testing menjadi satu data frame
FIPS22.train.xy$dataset <- "Training"
FIPS22.test.xy$dataset <- "Testing"
dim(FIPS22.test.xy)
dim(FIPS22.train.xy)
# Gabungkan kedua data
combined_data_GWRF <- rbind(FIPS22.train.xy, FIPS22.test.xy)

library(dplyr)
# Membuat kategori untuk kolom R-Square Lokal
combined_data_GWRF <- combined_data_GWRF %>%
  mutate(
    R2.Local.Cat = case_when(
      R2.GWRF <= 0.2 ~ "≤ 0.2",
      R2.GWRF > 0.2 & R2.GWRF <= 0.4 ~ "(0.2, 0.4]",
      R2.GWRF > 0.4 & R2.GWRF <= 0.6 ~ "(0.4, 0.6]",
      R2.GWRF > 0.6 & R2.GWRF <= 0.8 ~ "(0.6, 0.8]",
      R2.GWRF > 0.8 ~ "> 0.8"
    )
  ) %>%
  # Convert the category column into a factor for efficiency
  mutate(R2.Local.Cat = factor(R2.Local.Cat, levels = c("≤ 0.2", "(0.2, 0.4]", "(0.4, 0.6]", "(0.6, 0.8]", "> 0.8")))

# Display proportions as percentages
round(prop.table(table(combined_data_GWRF$R2.Local.Cat)) * 100,digits=3)

# Define colors using RColorBrewer palette
category_colors <- brewer.pal(n = 5, name = "Set3")  # Using Set3 palette with 5 colors

# Plot
r2.gwrf22 <- ggplot(data = combined_data_GWRF) +
  geom_sf(aes(fill = R2.Local.Cat, color = dataset), show.legend = TRUE) +
  scale_fill_manual(values = category_colors, name = NULL) + # Warna berdasarkan kategori
  scale_color_manual(values = c("Testing" = "red"), name = NULL) +
  labs(title = "R² Lokal GWRF") +
  theme_void() + # Menghilangkan grid dan background
  theme(
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
  ) +
  guides(
    fill = guide_legend(
      nrow = 1, byrow = TRUE,
      override.aes = list(color = category_colors) # Warna garis sesuai kategori
    ),
    color = guide_legend(
      override.aes = list(fill = "white", color = "red", size = 1.2) # Testing fill putih dengan garis merah
    )
  ) +
  coord_sf(xlim = c(105.5, 114.5), ylim = c(-8.8, -5.5)) # Zoom wilayah Jawa Timur

ggsave("D:/Kuliah S2/Thesis/Syntax/Eksplorasi/22_R2_GWRF1.png", plot = r2.gwrf22, dpi = 800,width = 10, height = 5,  units = "in")


#==== VISUALISASI DATA PREDIKSI
# Membuat kategori untuk kolom Pred.GWRF.Local
combined_data_GWRF <- combined_data_GWRF %>%
  mutate(
    Pred.GWRF.Local.Cat = case_when(
      Pred.GWRF.Local <= 3 ~ "< 3",
      Pred.GWRF.Local > 3 & Pred.GWRF.Local <= 5 ~ "3 - 5",
      Pred.GWRF.Local > 5 & Pred.GWRF.Local <= 7 ~ "5 - 7",
      Pred.GWRF.Local > 7 & Pred.GWRF.Local <= 10 ~ "7 - 10",
      Pred.GWRF.Local > 10 ~ "> 10"
    )
  ) %>%
  # Konversikan kolom kategori menjadi factor langsung untuk efisiensi
  mutate(Pred.GWRF.Local.Cat = factor(Pred.GWRF.Local.Cat, levels = c("< 3", "3 - 5", "5 - 7", "7 - 10", "> 10")))

# Definisikan warna kategori dengan cara yang lebih langsung
category_colors <- c(
  "< 3" = "#e0f4ff",  # Biru sangat muda
  "3 - 5" = "#c6dbef", # Biru muda
  "5 - 7" = "#6baed6", # Biru sedang
  "7 - 10" = "#2171b5", # Biru lebih gelap
  "> 10" = "#084594"   # Biru gelap lebih dalam
)

# Plot
pred.gwrf22 <- ggplot(data = combined_data_GWRF) +
  geom_sf(aes(fill = Pred.GWRF.Local.Cat, color = dataset), show.legend = TRUE) +
  scale_fill_manual(values = category_colors, name = NULL) + # Warna berdasarkan kategori
  scale_color_manual(values = c("Testing" = "red"), name =NULL) +
  labs(title = "Prediksi GWRF") +
  theme_void() + # Menghilangkan grid dan background
  theme(
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
  ) +
  guides(
    fill = guide_legend(
      nrow = 1, byrow = TRUE,
      override.aes = list(color = category_colors) # Warna garis sesuai kategori
    ),
    color = guide_legend(
      override.aes = list(fill = "white", color = "red", size = 1.2) # Testing fill putih dengan garis merah
    )
  ) +
  coord_sf(xlim = c(105.5, 114.5), ylim = c(-8.8, -5.5)) # Zoom wilayah Jawa Timur

# Simpan plot
ggsave("D:/Kuliah S2/Thesis/Syntax/Eksplorasi/22_Pred_GWRF1.png", plot = pred.gwrf22, dpi = 800,width = 10, height = 5,  units = "in")

#==== VISUALISASI DATA AKTUAL
combined_data_GWRF <- combined_data_GWRF %>%
  mutate(
    AKB_Cat = case_when(
      AKB <= 3 ~ "< 3",             # Untuk nilai < 3
      AKB > 3 & AKB <= 5 ~ "3 - 5",  # Untuk nilai antara 3 dan 5
      AKB > 5 & AKB <= 7 ~ "5 - 7",  # Untuk nilai antara 5 dan 7
      AKB > 7 & AKB <= 10 ~ "7 - 10",# Untuk nilai antara 7 dan 10
      AKB > 10 ~ "> 10"             # Untuk nilai lebih dari 10
    )
  ) %>%
  # Konversikan kolom kategori menjadi factor langsung untuk efisiensi
  mutate(AKB_Cat = factor(AKB_Cat, levels = c("< 3", "3 - 5", "5 - 7", "7 - 10", "> 10")))

# Plot
actual22 <- ggplot(data = combined_data_GWRF) +
  geom_sf(aes(fill = AKB_Cat, color = dataset), show.legend = TRUE) +
  scale_fill_manual(values = category_colors, name = NULL) + # Warna berdasarkan kategori
  scale_color_manual(values = c("Testing" = "red"), name =NULL) +
  labs(title = "Data Aktual") +
  theme_void() + # Menghilangkan grid dan background
  theme(
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
  ) +
  guides(
    fill = guide_legend(
      nrow = 1, byrow = TRUE,
      override.aes = list(color = category_colors) # Warna garis sesuai kategori
    ),
    color = guide_legend(
      override.aes = list(fill = "white", color = "red", size = 1.2) # Testing fill putih dengan garis merah
    )
  ) +
  coord_sf(xlim = c(105.5, 114.5), ylim = c(-8.8, -5.5)) # Zoom wilayah Jawa Timur

# Simpan plot
ggsave("D:/Kuliah S2/Thesis/Syntax/Eksplorasi/22_Data Aktual1.png", plot = actual22, dpi = 800,width = 10, height = 5,  units = "in")

#========= DISTRIBUSI VARIABEL IMPROTANCE =========# 
par(mfrow = c(1, 1))
#==== BAR PLOT
# Mengambil data dari l.VariableImportance
importance_df22 <- data.frame(
  Variable = rownames(model.gwrf22$LocalModelSummary$l.VariableImportance),
  Importance = model.gwrf22$LocalModelSummary$l.VariableImportance[, "Mean"]
) %>%
  arrange(desc(Importance)) %>%
  mutate(Rank = row_number())
importance_df22 

# Membuat grafik bar untuk local variable importance
barplot22 <- ggplot(importance_df22, aes(x = reorder(Variable, Importance), y = Importance)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(title = NULL, x = "Variabel", y = "%IncMSE") +
  theme_minimal() +
  coord_flip() + # Membalik sumbu x dan y untuk tampilan yang lebih baik
  theme(
    panel.grid = element_blank(),  # Menghapus grid
    axis.text.x = element_text(color = "black"),  # Mengatur warna teks pada sumbu X
    axis.text.y = element_text(color = "black")   # Mengatur warna teks pada sumbu Y
  )

ggsave("D:/Kuliah S2/Thesis/Syntax/Eksplorasi/22_GWRF_Barplot_Variabel_Importance.png", plot = barplot22, dpi = 800, width = 6, height = 4, units = "in")


#===== PLOT PIE : Peran variabel importance
# Gabungkan data WADMKK dengan nilai penting lokal dan tambahkan faktor utama
local_importance22 <- cbind(
  WADMKK = train22$WADMKK,
  model.gwrf22$Local.Variable.Importance,
  primary_factor = colnames(model.gwrf22$Local.Variable.Importance)[
    max.col(as.matrix(model.gwrf22$Local.Variable.Importance), ties.method = "first")
  ]
)
local_importance22[,c("WADMKK","primary_factor")]

# Hitung frekuensi dan proporsi faktor utama
factor_data22 <- as.data.frame(prop.table(table(local_importance22$primary_factor)) * 100) %>%
  rename(factor = Var1, proportion = Freq)
factor_data22

# Membuat pie chart dengan kolom yang benar
library(RColorBrewer)
pie(factor_data22$proportion,
    labels = paste0(factor_data22$factor, " (", round(factor_data22$proportion, 1), "%)"),
    main = "Local primary factor share",
    col = brewer.pal(n = length(factor_data22$factor), name = "Set2"),
    clockwise = TRUE)

#===== DISTRIBUSI SPATIAL : Local primary factor share
local_importance22 <- local_importance22 %>%
  mutate(x = train22$x, y = train22$y, geometry=train22$geometry)
local_importance22$primary_factor

# Data frame daerah untuk setiap kategori top3_vars
daerah_local_primary22 <- st_drop_geometry(local_importance22) %>%
  group_by(primary_factor) %>%
  summarize(
    regions = paste(WADMKK, collapse = ", ")  # Gabungkan nama daerah dalam setiap kategori
  ) %>%
  arrange(primary_factor)  # Mengurutkan berdasarkan kategori top3_vars
print(daerah_local_primary22, n = Inf)

#==== Visualisasi Plot Local Primary
gwrf.SPDF22<-train22[c(1:4,17)]
gwrf.SPDF22$primary_factor <- factor(local_importance22$primary_factor)
gwrf.sf <- st_as_sf(gwrf.SPDF22)
gwrf.SPDF22 <- as(gwrf.sf, "Spatial")

xlim_range <- c(105, 115)  # Tentukan nilai minimum dan maksimum untuk sumbu X
ylim_range <- c(-8.8, -5.5)

# Simpan dalam format png menggunakan perangkat grafis R
jpeg("D:/Kuliah S2/Thesis/Syntax/Eksplorasi/22_GWRF_Local Primary Factors1.png", 
     width = 10, height = 4, units = "in", res = 800)

# Plot menggunakan spplot
print(spplot(gwrf.SPDF22, "primary_factor", 
             main = "Spatial Distribution of Local Primary Factors", 
             sp.layout = list(state), 
             col = "transparent", 
             col.regions = rev(myPaletteRes(100)),
             auto.key = list(space = "bottom"), # Menempatkan legenda di bawah
             xlim = xlim_range,  # Menambahkan batas untuk sumbu X
             ylim = ylim_range)) 

# Tutup perangkat grafis
dev.off()


#======== DITRIBUSI SPASIAL PER VARIABEL
### Local feature importance (IncMSE)
gwrf.SPDF22$incMSE.PK <- model.gwrf22[["Local.Variable.Importance"]][["PK"]]
gwrf.SPDF22$incMSE.RLS <- model.gwrf22[["Local.Variable.Importance"]][["RLS"]]
gwrf.SPDF22$incMSE.PKKP <- model.gwrf22[["Local.Variable.Importance"]][["PKKP"]]
gwrf.SPDF22$incMSE.KP <- model.gwrf22[["Local.Variable.Importance"]][["KP"]]
gwrf.SPDF22$incMSE.SL <- model.gwrf22[["Local.Variable.Importance"]][["SL"]]
gwrf.SPDF22$incMSE.AML <- model.gwrf22[["Local.Variable.Importance"]][["AML"]]
gwrf.SPDF22$incMSE.IDL <- model.gwrf22[["Local.Variable.Importance"]][["IDL"]]
gwrf.SPDF22$incMSE.ASI <- model.gwrf22[["Local.Variable.Importance"]][["ASI"]]
gwrf.SPDF22$incMSE.BBLR <- model.gwrf22[["Local.Variable.Importance"]][["BBLR"]]
gwrf.SPDF22$incMSE.Rdokter <- model.gwrf22[["Local.Variable.Importance"]][["Rdokter"]]
gwrf.SPDF22$incMSE.Rpuskesmas <- model.gwrf22[["Local.Variable.Importance"]][["Rpuskesmas"]]

#gwrf.SPDF22@data[,c(1,12, 16, 17, 19)]

library(writexl)
write_xlsx(gwrf.SPDF22@data, path = "Variabel importance end.xlsx")
write_xlsx(combined_data_GWRF, path = "Resume All Model.xlsx")
#=== PLOTING VARIABLE IMPORTANCE LOCAL
library(RColorBrewer)
par(mfrow = c(1, 1))
col.palette.t <- colorRampPalette(c("#4C80A1", "#BEBADA","#8DD3C7", "#FFFFB3","#FB8072" ), 
                                  space="rgb", interpolate="linear")
# Daftar nama kolom yang akan dipetakan
columns <- c("incMSE.PK", "incMSE.RLS", "incMSE.PKKP", "incMSE.KP", 
             "incMSE.SL", "incMSE.AML", "incMSE.ASI", 
             "incMSE.Rdokter", "incMSE.Rpuskesmas", "incMSE.BBLR", "incMSE.IDL")

# Daftar judul untuk masing-masing kolom
titles <- c("Penduduk Miskin", "Rata-Rata Lama Sekolah", 
            "Prevelensi Ketidakcukupan Konsumsi Pangan", 
            "Ketahanan Pangan", "Akses Terhadap Layanan Sanitasi Layak", 
            "Akses Terhadap Sumber Air Minum Layak", "Asi Eksklusif", 
            "Rasio Dokter", "Rasio Puskesmas", "BBLR", 
            "Imunisasi Dasar Lengkap")

# Warna untuk kolom
col.palette.t <- colorRampPalette(c("#4C80A1", "#BEBADA","#8DD3C7", "#FFFFB3","#FB8072" ), 
                                  space="rgb", interpolate="linear")

# Iterasi untuk membuat dan menyimpan plot
for(i in 1:length(columns)) {
  # Buat nama file untuk setiap plot
  file_name <- paste0("D:/Kuliah S2/Thesis/Syntax/Eksplorasi/", i, "_", gsub(" ", "_", titles[i]), "22.png")
  
  # Simpan plot dalam format png
  jpeg(file_name, width = 10, height = 4, units = "in", res = 800)
  
  # Buat plot menggunakan spplot
  print(spplot(gwrf.SPDF22, columns[i], main = titles[i], 
               sp.layout = list(state),
               col = "transparent", 
               col.regions = rev(col.palette.t(100)), 
               xlim = xlim_range,  # Menambahkan batas untuk sumbu X
               ylim = ylim_range))
  
  # Tutup perangkat grafis
  dev.off()
}

#===== Local goodness of fit
#======= R-SQUARE
gwrf.SPDF22$Rsqrt <- model.gwrf22[["LGofFit"]][["LM_Rsq100"]]
gwrf.SPDF22$Rsqrt 
summary(gwrf.SPDF22$Rsqrt)

myPaletteRes <- colorRampPalette(c("lightseagreen","lightsteelblue1", "moccasin","hotpink", "red"))

# Simpan dalam format png menggunakan perangkat grafis R
jpeg("D:/Kuliah S2/Thesis/Syntax/Eksplorasi/22_GWRF_R2 Local OOB.jpg", 
     width = 10, height = 4, units = "in", res = 800)

# Plot menggunakan spplot
print(spplot(gwrf.SPDF22, "Rsqrt", main = "Local R2 (%)", 
             sp.layout = list(state),
             col = "transparent", 
             col.regions = rev(myPaletteRes(100)), 
             xlim = xlim_range,  # Menambahkan batas untuk sumbu X
             ylim = ylim_range))

# Tutup perangkat grafis
dev.off()

#====== RESIDUAL STANDARIZED
gwrf.SPDF22$Res <- model.gwrf22[["LGofFit"]]$LM_ResPred;gwrf.SPDF22$Res 
mean_res <- mean(gwrf.SPDF22$Res)
sd_res <- sd(gwrf.SPDF22$Res)

# Tambahkan standardized residuals ke data
gwrf.SPDF22$Standardized_Res <- (gwrf.SPDF22$Res - mean_res) / sd_res;gwrf.SPDF22$Standardized_Res

# Simpan dalam format png menggunakan perangkat grafis R
jpeg("D:/Kuliah S2/Thesis/Syntax/Eksplorasi/22_GWRF_Standardized_Residual.png", 
     width = 10, height = 4, units = "in", res = 800)

# Plot menggunakan spplot
print(spplot(gwrf.SPDF22, "Standardized_Res", main = "Local Standardized Residual", 
             sp.layout = list(state),
             col = "transparent", 
             col.regions = rev(myPaletteRes (100)), 
             xlim = xlim_range,  # Menambahkan batas untuk sumbu X
             ylim = ylim_range))

# Tutup perangkat grafis
dev.off()


#====== LOCAL MORAN'S I ======
library(spdep)
coords22 <- coordinates(gwrf.SPDF22)
neighbors22 <- knearneigh(coords22, k = 5) # Tetapkan minimal 5 tetangga
weights22 <- nb2listw(knn2nb(neighbors22), style = "W")
# Uji Moran's I
moran_result1 <- moran.test(gwrf.SPDF22$Standardized_Res , weights22)
print(moran_result1)

# 2. Hitung Local Moran's I
moran_local <- localmoran(gwrf.SPDF22$Standardized_Res , weights22)


# 3. Tambahkan hasil Local Moran's I ke objek spasial
gwrf.SPDF22$Local_Moran_I <- moran_local[, 1]  # Indeks Moran lokal
gwrf.SPDF22$p_value <- moran_local[, 5]       # Nilai p dari tes signifikansi

# 4. Visualisasi peta Local Moran's I
# Membuat peta nilai Local Moran's I
jpeg("D:/Kuliah S2/Thesis/Syntax/Eksplorasi/22_GWRF_Local_Moran Indeks.png", 
     width = 10, height = 4, units = "in", res = 800)

print(spplot(gwrf.SPDF22, "Local_Moran_I", main = "Local Moran's I of Residuals",
             sp.layout = list(state),
             col = "transparent", 
             col.regions = rev(myPaletteRes(100)), 
             xlim = xlim_range,  # Menambahkan batas untuk sumbu X
             ylim = ylim_range))

# Tutup perangkat grafis
dev.off()

# 5. (Opsional) Visualisasi peta p-value
jpeg("D:/Kuliah S2/Thesis/Syntax/Eksplorasi/22_GWRF_Local_Moran_pvalue.png", 
     width = 10, height = 4, units = "in", res = 800)

print(spplot(gwrf.SPDF22, "p_value", main = "P-Value of Local Moran's I",
             sp.layout = list(state),
             col = "transparent", 
             col.regions = rev(myPaletteRes(100)), 
             xlim = xlim_range,  # Menambahkan batas untuk sumbu X
             ylim = ylim_range))

dev.off()

#==== Cluster Moran
# Definisikan threshold untuk menentukan High dan Low (z-score 0 sebagai median)
gwrf.SPDF22$Cluster <- with(gwrf.SPDF22, 
                            ifelse(gwrf.SPDF22$p_value < 0.05 & gwrf.SPDF22$Res > 0 & gwrf.SPDF22$Local_Moran_I > 0, "High-High",
                                   ifelse(gwrf.SPDF22$p_value < 0.05 & gwrf.SPDF22$Res > 0 & gwrf.SPDF22$Local_Moran_I < 0, "High-Low",
                                          ifelse(gwrf.SPDF22$p_value< 0.05 & gwrf.SPDF22$Res < 0 & gwrf.SPDF22$Local_Moran_I > 0, "Low-Low",
                                                 ifelse(gwrf.SPDF22$p_value < 0.05 & gwrf.SPDF22$Res < 0 & gwrf.SPDF22$Local_Moran_I < 0, "Low-High", "Not Significant")))))


# Konversi ke faktor untuk mempermudah visualisasi
gwrf.SPDF22$Cluster <- factor(gwrf.SPDF22$Cluster, levels = c("High-High", "High-Low", "Low-Low", "Low-High",  "Not Significant"))

# Menghitung frekuensi faktor utama
factor_counts1 <- table(gwrf.SPDF22@data$Cluster)
factor_counts1

library(RColorBrewer)
# Pilih warna untuk kategori cluster
cluster_colors <- c("red","pink","lightblue", "blue", "gray") # HH, HL, LH, LL, None

jpeg("D:/Kuliah S2/Thesis/Syntax/Eksplorasi/22_GWRF_Local_Moran_Clusters1.png", 
     width = 10, height = 3.5, units = "in", res = 800)

print(spplot(gwrf.SPDF22, "Cluster",
             sp.layout = list(state),
             col.regions = cluster_colors,
             main = "Spatial Clusters (Local Moran's I)", 
             xlim = xlim_range,  # Menambahkan batas untuk sumbu X
             ylim = ylim_range))

dev.off()


#======================= RANDOM FOREST
library(randomForest)
# Definisikan target dan predictors
target <- "AKB"
predictors <- c("PK", "RLS", "PKKP", "KP", "SL", "AML", 
                "IDL", "ASI", "BBLR", "Rdokter", "Rpuskesmas")

# Buat formula
formula_rf <- as.formula(paste(target, "~", paste(predictors, collapse = "+")))

model_rf <- randomForest(
  formula = formula_rf,
  data = train.df22,
  ntrees = 500,
  mtry = 4,  # Proporsi kolom
  nodesize = 5,                           # Minimum rows
  max.depth = 10,                              # Max depth
  importance = T,                 # Variable importance
  keep.forest = T,
  seed = 1234567
)

summary(model_rf)


############### POHON KEPUTUSAN RF
# Mengambil jumlah pohon dalam model RF  
num_trees <- model_rf$ntree  

# Variabel untuk menyimpan informasi tentang pohon paling sederhana  
simple_tree_info <- list(tree = NULL, depth = Inf, node_count = Inf, tree_index = NA)  

# Loop melalui setiap pohon untuk mencari yang paling sederhana  
for (i in 1:num_trees) {  
  tree <- getTree(model_rf, k = i, labelVar = TRUE)  
  
  # Hitung kedalaman pohon  
  depth <- max(tree$left.daughter[tree$status == -1], na.rm = TRUE)  
  node_count <- nrow(tree)  
  
  # Cek jika ini pohon dengan konstruk paling sederhana  
  if (depth < simple_tree_info$depth || (depth == simple_tree_info$depth && node_count < simple_tree_info$node_count)) {  
    simple_tree_info$tree <- tree  
    simple_tree_info$depth <- depth  
    simple_tree_info$node_count <- node_count  
    simple_tree_info$tree_index <- i  # Menyimpan index pohon yang paling sederhana  
  }  
}  

# Menampilkan informasi tentang pohon paling sederhana yang ditemukan  
if (!is.null(simple_tree_info$tree)) {  
  cat("Pohon keputusan yang paling sederhana adalah pohon ke:", simple_tree_info$tree_index, "\n")  
  cat("Kedalaman pohon:", simple_tree_info$depth, "\n")  
  cat("Jumlah node:", simple_tree_info$node_count, "\n")  
}

# Buat fungsi untuk mengonversi ke format DOT dan menggambar  
write_dot <- function(tree_data) {  
  dot_string <- "digraph Tree {\n"  
  dot_string <- paste(dot_string, "  node [shape = rectangle, style = filled, fillcolor = lightblue]\n")  
  
  for (i in seq_len(nrow(tree_data))) {  
    label <- ifelse(tree_data$status[i] == -1,  
                    paste( round(as.numeric(tree_data$prediction[i]), 2)),  
                    paste(tree_data$`split var`[i], "<", round(as.numeric(tree_data$`split point`[i]), 2)))  
    
    dot_string <- paste(dot_string, sprintf("  %d [label=\"%s\"];\n", i, label))  
    
    if (tree_data$`left daughter`[i] != 0) {  
      dot_string <- paste(dot_string, sprintf("  %d -> %d [label = 'Yes'];\n", i, tree_data$`left daughter`[i]))  
    }  
    if (tree_data$`right daughter`[i] != 0) {  
      dot_string <- paste(dot_string, sprintf("  %d -> %d [label = 'No'];\n", i, tree_data$`right daughter`[i]))  
    }  
  }  
  
  dot_string <- paste(dot_string, "}\n")  
  return(dot_string)  
}  

# Generate the DOT string from the most simple tree  
if (!is.null(simple_tree_info$tree)) {  
  dot_data <- write_dot(simple_tree_info$tree)  
  
  # Visualize with grViz  
  library(DiagrammeR)  
  grViz(dot_data)  
} else {  
  cat("Tidak ada pohon yang ditemukan.")  
} 

# Ekstrak pohon pertama (tree ke-1)
tree1 <- getTree(model_rf, k = 251, labelVar = TRUE)

write_dot <- function(tree_data) {  
  # Membuat grafik DOT  
  dot_string <- "digraph Tree {\n"  
  dot_string <- paste(dot_string, "  node [shape = rectangle, style = filled, fillcolor = lightblue]\n")  
  
  for (i in seq_len(nrow(tree_data))) {  
    # Mendapatkan label  
    label <- ifelse(tree_data$status[i] == -1,  
                    paste("Prediction: ", round(as.numeric(tree_data$prediction[i]), 2)),  
                    paste(tree_data$`split var`[i], "<", round(as.numeric(tree_data$`split point`[i]), 2)))  
    
    # Menambahkan node ke STRING  
    dot_string <- paste(dot_string, sprintf("  %d [label=\"%s\"];\n", i, label))  
    
    # Menambahkan edges  
    if (tree_data$`left daughter`[i] != 0) {  
      dot_string <- paste(dot_string, sprintf("  %d -> %d [label = 'Yes'];\n", i, tree_data$`left daughter`[i]))  
    }  
    
    if (tree_data$`right daughter`[i] != 0) {  
      dot_string <- paste(dot_string, sprintf("  %d -> %d [label = 'No'];\n", i, tree_data$`right daughter`[i]))  
    }  
  }  
  
  dot_string <- paste(dot_string, "}\n")  
  return(dot_string)  
}  

dot_data <- write_dot(tree1)  

# Visualize with grViz  
grViz(dot_data)


#============== VARIABEL IMPORTANCE
var.imp <- importance(model_rf)

importance_values <- importance(model_rf, type = 1)  # type = 1 untuk IncMSE
importanceRF_df <- data.frame(
  Variable = rownames(importance_values),
  IncMSE = importance_values[, "%IncMSE"]
)
importanceRF_df <- importanceRF_df[order(importanceRF_df$IncMSE, decreasing = TRUE), ]

var.impRF <-ggplot(importanceRF_df, aes(x = reorder(Variable, IncMSE), y = IncMSE)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +  # Membalik sumbu untuk horizontal bar
  labs(
    title = NULL,
    x = "Variabel",
    y = "%IncMSE"
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),  # Menghapus grid
    axis.text.x = element_text(color = "black"),  # Mengatur warna teks pada sumbu X
    axis.text.y = element_text(color = "black")   # Mengatur warna teks pada sumbu Y
  )

ggsave("D:/Kuliah S2/Thesis/Syntax/Eksplorasi/22_RF_Variabel Importance.jpg", plot = var.impRF, dpi = 800,width = 6, height = 4,  units = "in")

#==== PREDIKSI
#dari model RF
FIPS22.train.xy$Pred.RF.Global1 <- model_rf$predicted
pred_rf<- predict(model_rf ,test.df22)
FIPS22.test.xy$Pred.RF.Global1 <- pred_rf

# Nilai target aktual dari data pelatihan
y_actual <- train.df22$AKB
y_pred <- predict(model_rf, train.df22)
y_mean <- mean(y_actual)
TSS <- sum((y_actual - y_mean)^2)
RSS <- sum((y_actual - y_pred)^2)
R_squared_not_oob <- 1 - (RSS / TSS)
R_squared_not_oob

cat('RF RMSE:', round(sqrt(mean((FIPS22.train.xy$Pred.RF.Global1-FIPS22.train.xy$AKB)^2 , na.rm = TRUE)), digits=3), '\n')
cat('RF MAE:', round(mean(abs(FIPS22.train.xy$Pred.RF.Global1-FIPS22.train.xy$AKB) , na.rm = TRUE ), digits=3), '\n')
cat('RF R2:', round(R_squared_test, digits=3), '\n')

# Nilai target aktual dari data testing
y_actual_test <- test.df22$AKB
y_pred_test <- predict(model_rf, test.df22)
y_mean_test <- mean(y_actual_test)
TSS_test <- sum((y_actual_test - y_mean_test)^2)
RSS_test <- sum((y_actual_test - y_pred_test)^2)
R_squared_test <- 1 - (RSS_test / TSS_test)
R_squared_test

cat('RF Test RMSE:', round(sqrt(mean((FIPS22.test.xy$Pred.RF.Global1-FIPS22.test.xy$AKB)^2 , na.rm = TRUE)), digits=3), '\n')
cat('RF Test MAE:', round(mean(abs(FIPS22.test.xy$Pred.RF.Global1-FIPS22.test.xy$AKB) , na.rm = TRUE ), digits=3), '\n')
cat('RF Test R2:', round(R_squared_test, digits=3), '\n')
round(summary(lm(AKB~Pred.RF.Global1,FIPS22.test.xy))$r.squared, digits=3)

# Hitung nilai R2 dan RMSE untuk dimasukan ke dalam plot Testing
rf_rmse <- round(sqrt(mean((FIPS22.test.xy$Pred.RF.Global1 - FIPS22.test.xy$AKB)^2, na.rm = TRUE)), digits = 3)
rf_r2 <- round(R_squared_test, digits=3)
rf_mae <- round(mean(abs(FIPS22.test.xy$Pred.RF.Global1-FIPS22.test.xy$AKB)), digits = 3)

#=== Ploting actual vs prediction Testing
p.test.RF22 <- ggplot(data = FIPS22.test.xy, aes(x = AKB, y = Pred.RF.Global1)) + 
  geom_point(size = 2.0, color = "blue") +  # Titik berwarna biru
  geom_smooth(method = "lm", se = FALSE, colour = "black") +  # Garis regresi berwarna hitam
  labs(x = "Observed", y = "Predicted", 
       title = "Observed vs Predicted RF") +
  geom_label(aes(
    x = min(FIPS22.test.xy$AKB, na.rm = TRUE) + 0.1 * diff(range(FIPS22.test.xy$AKB, na.rm = TRUE)), 
    y = max(FIPS22.test.xy$Pred.RF.Global1, na.rm = TRUE) - 0.1 * diff(range(FIPS22.test.xy$Pred.RF.Global1, na.rm = TRUE)), 
    label = paste0("MAE: ", rf_mae, "\nRMSE: ", rf_rmse, "\nR²: ", rf_r2)),
    hjust = 0.5, vjust = 0.2, size = 4.5, fill = "white", color = "black", label.size = 0.1) +  # Warna kotak putih
  theme_minimal() +  # Latar belakang putih dengan tema minimal
  theme(legend.position = "top")  # Legend di atas

# Menampilkan plot
print(p.test.RF22)

ggsave("D:/Kuliah S2/Thesis/Syntax/Eksplorasi/22_RF_Residual Eval Testing.png", plot = p.test.RF22, dpi = 800,width = 5.2, height = 5.2,  units = "in")


#=== Ploting actual vs prediction Training
trainrf_rmse <- round(sqrt(mean((FIPS22.train.xy$Pred.RF.Global1 - FIPS22.train.xy$AKB)^2, na.rm = TRUE)), digits = 3)
trainrf_r2 <- round(R_squared_not_oob, digits = 3)
trainrf_mae <- round(mean(abs(FIPS22.train.xy$Pred.RF.Global1-FIPS22.train.xy$AKB)), digits = 3)

p.test.RF22train <- ggplot(data = FIPS22.train.xy, aes(x = AKB, y = Pred.RF.Global1)) + 
  geom_point(size = 2.0, color = "blue") +  # Titik berwarna biru
  geom_smooth(method = "lm", se = FALSE, colour = "black") +  # Garis regresi berwarna hitam
  labs(x = "Observed", y = "Predicted", 
       title = "Observed vs Predicted RF") +
  geom_label(aes(
    x = min(FIPS22.train.xy$AKB, na.rm = TRUE) + 0.1 * diff(range(FIPS22.train.xy$AKB, na.rm = TRUE)), 
    y = max(FIPS22.train.xy$Pred.RF.Global1, na.rm = TRUE) - 0.1 * diff(range(FIPS22.train.xy$Pred.RF.Global1, na.rm = TRUE)), 
    label = paste0("MAE: ", trainrf_mae, "\nRMSE: ", trainrf_rmse, "\nR²: ", trainrf_r2)),
    hjust = 0.5, vjust = 0.2, size = 4.5, fill = "white", color = "black", label.size = 0.1) +  # Warna kotak putih
  theme_minimal() +  # Latar belakang putih dengan tema minimal
  theme(legend.position = "top")  # Legend di atas

# Menampilkan plot
print(p.test.RF22train)

ggsave("D:/Kuliah S2/Thesis/Syntax/Eksplorasi/22_RF_Residual Eval Training.png", plot = p.test.RF22train, dpi = 800,width = 5.2, height = 5.2,  units = "in")


#===== DISTRIBUSI SPASIAL PREDIKSI FULL
# Gabungkan data training dan testing menjadi satu data frame
FIPS22.train.xy$dataset <- "Training"
FIPS22.test.xy$dataset <- "Testing"
combined_data_GWRF
# Gabungkan kedua data
combined_data_GWRF <- rbind(FIPS22.train.xy, FIPS22.test.xy)
summary(combined_data$Pred.RF.Global1)
# Membuat kategori untuk kolom Pred.RF.Global1
combined_data_GWRF <- combined_data_GWRF %>%
  mutate(
    Pred.RF.Global1.Cat = case_when(
      Pred.RF.Global1 <= 3 ~ "< 3",
      Pred.RF.Global1 > 3 & Pred.RF.Global1 <= 5 ~ "3 - 5",
      Pred.RF.Global1 > 5 & Pred.RF.Global1 <= 7 ~ "5 - 7",
      Pred.RF.Global1 > 7 & Pred.RF.Global1 <= 10 ~ "7 - 10",
      Pred.RF.Global1 > 10 ~ "> 10"
    )
  ) %>%
  # Konversikan kolom kategori menjadi factor langsung untuk efisiensi
  mutate(Pred.RF.Global1.Cat = factor(Pred.RF.Global1.Cat, levels = c("< 3", "3 - 5", "5 - 7", "7 - 10", "> 10")))

# Definisikan warna kategori dengan cara yang lebih langsung
category_colors <- c(
  "3 - 5" = "#c6dbef", # Biru muda
  "5 - 7" = "#6baed6", # Biru sedang
  "7 - 10" = "#2171b5", # Biru lebih gelap
  "> 10" = "#084594"   # Biru gelap lebih dalam
)

# Plot
pred.rf22 <- ggplot(data = combined_data_GWRF) +
  geom_sf(aes(fill = Pred.RF.Global1.Cat, color = dataset), show.legend = TRUE) +
  scale_fill_manual(values = category_colors, name = NULL) + # Warna berdasarkan kategori
  scale_color_manual(values = c("Testing" = "red"), name =NULL) +
  labs(title = "Prediksi RF") +
  theme_void() + # Menghilangkan grid dan background
  theme(
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
  ) +
  guides(
    fill = guide_legend(
      nrow = 1, byrow = TRUE,
      override.aes = list(color = category_colors) # Warna garis sesuai kategori
    ),
    color = guide_legend(
      override.aes = list(fill = "white", color = "red", size = 1.2) # Testing fill putih dengan garis merah
    )
  ) +
  coord_sf(xlim = c(105.5, 114.5), ylim = c(-8.8, -5.5)) # Zoom wilayah Jawa Timur


# Simpan plot
ggsave("D:/Kuliah S2/Thesis/Syntax/Eksplorasi/22_RF_Pred1.png", plot = pred.rf22, dpi = 800,width = 10, height = 5,  units = "in")


#============= PDP
predictors <- c("PK", "RLS", "PKKP", "KP", "SL", "AML", 
                "IDL", "ASI", "BBLR", "Rdokter", "Rpuskesmas")

jpeg(
  filename = "D:/Kuliah S2/Thesis/Syntax/Eksplorasi/22_RF_PDP new.png",
  width = 18, height = 15, units = "in", res = 800
)

# Set up a 4x3 grid of plots
par(mfrow = c(4, 3))

# Create partial plots with larger text
partialPlot(model_rf, pred.data = train.df22, x.var = "BBLR", 
            n.pt = min(length(unique(train.df22$BBLR)), 51),
            rug = TRUE, 
            xlab = "BBLR", ylab = "Predicted AKB", 
            main = "(a) %BBLR", 
            cex.axis = 1.5, cex.lab = 1.5, cex.main = 2)

partialPlot(model_rf, pred.data = train.df22, x.var = "SL", 
            n.pt = min(length(unique(train.df22$SL)), 51),
            rug = TRUE, 
            xlab = "SL", ylab = "Predicted AKB", 
            main = "(b) %Akses Terhadap Layanan Sanitasi Layak",
            cex.axis = 1.5, cex.lab = 1.5, cex.main = 2)

partialPlot(model_rf, pred.data = train.df22, x.var = "AML", 
            n.pt = min(length(unique(train.df22$AML)), 51),
            rug = TRUE, 
            xlab = "AML", ylab = "Predicted AKB", 
            main = "(c) %Akses Terhadap Sumber Air Minum Layak",
            cex.axis = 1.5, cex.lab = 1.5, cex.main = 2)

partialPlot(model_rf, pred.data = train.df22, x.var = "RLS", 
            n.pt = min(length(unique(train.df22$RLS)), 51),
            rug = TRUE, 
            xlab = "RLS", ylab = "Predicted AKB", 
            main = "(d) Rata-rata Lama Sekolah",
            cex.axis = 1.5, cex.lab = 1.5, cex.main = 2)

partialPlot(model_rf, pred.data = train.df22, x.var = "PK", 
            n.pt = min(length(unique(train.df22$PK)), 51),
            rug = TRUE, 
            xlab = "PK", ylab = "Predicted AKB", 
            main = "(e) %Penduduk Miskin",
            cex.axis = 1.5, cex.lab = 1.5, cex.main = 2)

partialPlot(model_rf, pred.data = train.df22, x.var = "KP", 
            n.pt = min(length(unique(train.df22$KP)), 51),
            rug = TRUE, 
            xlab = "KP", ylab = "Predicted AKB", 
            main = "(f) Ketahanan Pangan",
            cex.axis = 1.5, cex.lab = 1.5, cex.main = 2)

partialPlot(model_rf, pred.data = train.df22, x.var = "PKKP", 
            n.pt = min(length(unique(train.df22$PKKP)), 51),
            rug = TRUE, 
            xlab = "PKKP", ylab = "Predicted AKB", 
            main = "(g) Prevalensi Ketidakcukupan Konsumsi Pangan",
            cex.axis = 1.5, cex.lab = 1.5, cex.main = 1.5)

partialPlot(model_rf, pred.data = train.df22, x.var = "ASI", 
            n.pt = min(length(unique(train.df22$ASI)), 51),
            rug = TRUE, 
            xlab = "ASI", ylab = "Predicted AKB", 
            main = "(h) %ASI Eksklusif",
            cex.axis = 1.5, cex.lab = 1.5, cex.main = 2)

partialPlot(model_rf, pred.data = train.df22, x.var = "Rdokter", 
            n.pt = min(length(unique(train.df22$Rdokter)), 51),
            rug = TRUE, 
            xlab = "Rdokter", ylab = "Predicted AKB", 
            main = "(i) Rasio Dokter",
            cex.axis = 1.5, cex.lab = 1.5, cex.main = 2)

partialPlot(model_rf, pred.data = train.df22, x.var = "Rpuskesmas", 
            n.pt = min(length(unique(train.df22$Rpuskesmas)), 51),
            rug = TRUE, 
            xlab = "Rpuskesmas", ylab = "Predicted AKB", 
            main = "(j) Rasio Puskesmas",
            cex.axis = 1.5, cex.lab = 1.5, cex.main = 2)

partialPlot(model_rf, pred.data = train.df22, x.var = "IDL", 
            n.pt = min(length(unique(train.df22$IDL)), 51),
            rug = TRUE, 
            xlab = "IDL", ylab = "Predicted AKB", 
            main = "(k) %Imunisasi Dasar Lengkap",
            cex.axis = 1.5, cex.lab = 1.5, cex.main = 2)

# Close the graphical device
dev.off()

#====================== MODEL GWR
library(GWmodel)
coords <- st_coordinates(train22) 
spatial_train22 <- as_Spatial(train22)
spatial_test22 <- as_Spatial(test22) 
class(spatial_train22)

spatial_train22@data[,7:17]<- scale(spatial_train22@data[,7:17])
spatial_test22@data[,7:17]<- scale(spatial_test22@data[,7:17])
str(spatial_train22@data)

#==== BANDWIDTH ADATIF
# BANDWIDTH BISQUARE
bw.Abi22 <- bw.gwr(AKB ~   PK + RLS + KP+PKKP  + SL + AML +  IDL + ASI + BBLR + Rdokter + Rpuskesmas, 
                   data = spatial_train22, approach = "AICc", adaptive = T, kernel="bisquare")

# MODEL
gwr_model.Abi22 <- gwr.basic(AKB ~   PK + RLS + KP+PKKP  + SL + AML +  IDL + ASI + BBLR + Rdokter + Rpuskesmas , 
                             data = spatial_train22,
                             bw = bw.Abi22, adaptive = T,kernel="bisquare")
gwr_model.Abi22$SDF
eval.Abi22 <- data.frame(gwr_model.Abi22[["GW.diagnostic"]])
eval.Abi22$RMSE <-sqrt(gwr_model.Abi22[["GW.diagnostic"]][["RSS.gw"]]/nrow(data22))
eval.Abi22 

# BANDWIDTH GAUSSIAN
bw.Agau22 <- bw.gwr(AKB ~   PK + RLS + KP+PKKP  + SL + AML +  IDL + ASI + BBLR + Rdokter + Rpuskesmas , data = spatial_train22, 
                    approach = "AICc", adaptive = T, kernel="gaussian")
# MODEL
gwr_model.Agau22 <- gwr.basic(AKB ~   PK + RLS + KP+PKKP  + SL + AML +  IDL + ASI + BBLR + Rdokter + Rpuskesmas , data = spatial_train22,
                              bw = bw.Agau22, adaptive = T,kernel="gaussian")
eval.Agau22 <- data.frame(gwr_model.Agau22[["GW.diagnostic"]])
eval.Agau22$RMSE <-sqrt(gwr_model.Agau22[["GW.diagnostic"]][["RSS.gw"]]/nrow(data22))
eval.Agau22 

# BANDWIDTH TRICUBE
bw.Atri22 <- bw.gwr(AKB ~   PK + RLS + KP+PKKP  + SL + AML +  IDL + ASI + BBLR + Rdokter + Rpuskesmas , data = spatial_train22, 
                    approach = "AICc", adaptive = T, kernel="tricube")
# MODEL
gwr_model.Atri22 <- gwr.basic(AKB ~   PK + RLS +  KP+PKKP+ SL + AML +  IDL + ASI + BBLR + Rdokter + Rpuskesmas , data = spatial_train22,
                              bw = bw.Atri22, adaptive = T,kernel="tricube")
eval.Atri22 <- data.frame(gwr_model.Atri22[["GW.diagnostic"]])
eval.Atri22$RMSE <-sqrt(gwr_model.Atri22[["GW.diagnostic"]][["RSS.gw"]]/nrow(data22))
eval.Atri22 

# BANDWIDTH EXPONENTIAL
bw.Aexp22 <- bw.gwr(AKB ~   PK + RLS +  PKKP + SL + AML +  IDL + ASI + BBLR + Rdokter + Rpuskesmas , data = spatial_train22, 
                    approach = "AICc", adaptive = T, kernel="exponential")
# MODEL
gwr_model.Aexp22 <- gwr.basic(AKB ~   PK + RLS +  PKKP + SL + AML +  IDL + ASI + BBLR + Rdokter + Rpuskesmas , data = spatial_train22,
                              bw = bw.Aexp22, adaptive = T,kernel="exponential")
eval.Aexp22 <- data.frame(gwr_model.Aexp22[["GW.diagnostic"]])
eval.Aexp22$RMSE <-sqrt(gwr_model.Aexp22[["GW.diagnostic"]][["RSS.gw"]]/nrow(data22))
eval.Aexp22 

kernelmod <- c("bisquare","gaussian","tricube","exponential")
hasil.gwr22 <- rbind(eval.Abi22, eval.Agau22, eval.Atri22, eval.Aexp22)
hasil.gwr22$kernel <- kernelmod
hasil.gwr22 <- hasil.gwr22[, c("kernel", setdiff(names(hasil.gwr22), "kernel"))]
hasil.gwr22

best_model22 <- hasil.gwr22[which.min(hasil.gwr22$RMSE), ] # Prioritaskan AICc
best_model22

best_model22 <- hasil.gwr22[which.max(hasil.gwr22$gw.R2), ] # Prioritaskan AICc
best_model22

best_kernel22 <- best_model22[,1]
best_kernel22
#========================================== MODEL TERBAIK 
best.bw22 <- bw.gwr(AKB ~   PK + RLS + KP+PKKP  + SL + AML +  IDL + ASI + BBLR + Rdokter + Rpuskesmas  , data = spatial_train22, 
                    approach = "AICc", adaptive = T, kernel=best_kernel22)
# MODEL
best.gwr22<- gwr.basic(AKB ~   PK + RLS + KP+PKKP  + SL + AML +  IDL + ASI + BBLR + Rdokter + Rpuskesmas  , data = spatial_train22,
                       bw = best.bw22, adaptive = T,kernel=best_kernel22)
best.gwr22

#########
#================= VISUALISASI
# Menambahkan koefisien untuk setiap variabel ke dalam data_sf
data_gwr <- spatial_train22[,c(1:5)]
data_gwr$PK_coef <- best.gwr22$SDF$PK           # Koefisien untuk PK
data_gwr$RLS_coef <- best.gwr22$SDF$RLS         # Koefisien untuk RLS
data_gwr$KP_coef <- best.gwr22$SDF$KP           # Koefisien untuk KP
data_gwr$PKKP_coef <- best.gwr22$SDF$PKKP       # Koefisien untuk PKKP
data_gwr$SL_coef <- best.gwr22$SDF$SL           # Koefisien untuk Sanitasi
data_gwr$AML_coef <- best.gwr22$SDF$AML         # Koefisien untuk AML
data_gwr$IDL_coef <- best.gwr22$SDF$IDL         # Koefisien untuk IDL
data_gwr$ASI_coef <- best.gwr22$SDF$ASI         # Koefisien untuk ASI
data_gwr$BBLR_coef <- best.gwr22$SDF$BBLR       # Koefisien untuk BBLR
data_gwr$Rdokter_coef <- best.gwr22$SDF$Rdokter # Koefisien untuk Rdokter
data_gwr$Rpuskesmas_coef <- best.gwr22$SDF$Rpuskesmas # Koefisien untuk Rpuskesmas
data_gwr$RSq <- best.gwr22$SDF$Local_R2

p.value <- gwr.t.adjust(best.gwr22)$results$p
head(p.value)
data_gwr <- cbind(data_gwr,p.value)
data_gwr$Standres <- (data_gwr$res - mean(data_gwr$res)) / sd(data_gwr$res)

data_gwr@data

local_main<- cbind(
  WADMKK = data_gwr$WADMKK,
  data_gwr@data[,c(6:16)],
  main_factor = colnames(data_gwr@data[,c(6:16)])[
    max.col(as.matrix(data_gwr@data[,c(6:16)]), ties.method = "first")
  ]
)

# Hitung frekuensi dan proporsi faktor utama
mainfactor_gwr <- as.data.frame(prop.table(table(local_main$main_factor)) * 100) %>%
  rename(factor = Var1, proportion = Freq)
mainfactor_gwr

# Mengakses slot @data untuk mengubah data
data_gwr$local_primary <- as.factor(local_main[,c("main_factor")])
data_gwr@data$local_primary <- gsub("_coef$", "", data_gwr@data$local_primary)
data_gwr@data$local_primary <- as.factor(data_gwr@data$local_primary)

# Mengecek hasil
head(data_gwr@data)

#===== FAKTOR PRIMER LOKAL GWR
jpeg("D:/Kuliah S2/Thesis/Syntax/Eksplorasi/22_GWR_Local Primary Factors1.png", 
     width = 10, height = 4, units = "in", res = 800)
print(spplot(data_gwr, "local_primary", 
             main = NULL, 
             sp.layout = list(state), 
             col = "transparent", 
             col.regions = rev(myPaletteRes(100)),
             auto.key = list(space = "bottom", 
                             title = NULL,
                             columns = 2, 
                             labels = legend_names),  # Menetapkan label baru
             xlim = xlim_range,  # Batas sumbu X
             ylim = ylim_range)) # Batas sumbu Y
dev.off()

#=== variabel BBLR
jpeg("D:/Kuliah S2/Thesis/Syntax/Eksplorasi/22_GWR_BBLR.png", 
     width = 10, height = 4, units = "in", res = 800)
print(spplot(data_gwr, "BBLR_coef", 
             main = "Distribusi Spasial BBLR", 
             sp.layout = list(state), 
             col = "transparent", 
             col.regions = rev(col.palette.t(100)),
             auto.key = list(space = "bottom"),
             xlim = xlim_range,  # Menambahkan batas untuk sumbu X
             ylim = ylim_range)) # Menempatkan legenda di bawah
dev.off()

#=== variabel puskesmas
# Simpan dalam format png menggunakan perangkat grafis R
jpeg("D:/Kuliah S2/Thesis/Syntax/Eksplorasi/22_GWR_Puskesmas.png", 
     width = 10, height = 4, units = "in", res = 800)

# Plot menggunakan spplot
print(spplot(data_gwr, "Rpuskesmas_coef", 
             main = "Distribusi Spasial Rasio Puskesmas", 
             sp.layout = list(state), 
             col = "transparent", 
             col.regions = rev(col.palette.t(100)),
             auto.key = list(space = "bottom"),
             xlim = xlim_range,  # Menambahkan batas untuk sumbu X
             ylim = ylim_range)) # Menempatkan legenda di bawah

# Tutup perangkat grafis
dev.off()

#=== variabel Rsquare
xlim_range <- c(105, 115)  # Tentukan nilai minimum dan maksimum untuk sumbu X
ylim_range <- c(-8.8, -5.5)
# Simpan dalam format png menggunakan perangkat grafis R
jpeg("D:/Kuliah S2/Thesis/Syntax/Eksplorasi/22_GWR_Rsquare.png", 
     width = 10, height = 4, units = "in", res = 800)

# Plot menggunakan spplot
print(spplot(data_gwr, "RSq", 
             main = "Distribusi Spasial R-square", 
             sp.layout = list(state), 
             col = "transparent", 
             col.regions = rev(col.palette.t(100)),
             auto.key = list(space = "bottom"),
             xlim = xlim_range,  # Menambahkan batas untuk sumbu X
             ylim = ylim_range))# Menempatkan legenda di bawah

# Tutup perangkat grafis
dev.off()


######=========== PREDIKSI
library(GWmodel)

# Prediksi untuk data testing
best.gwr22<- gwr.basic(AKB ~   PK + RLS + KP+PKKP  + SL + AML +  IDL + ASI + BBLR + Rdokter + Rpuskesmas  , data = spatial_train22,
                       bw = best.bw22, adaptive = T,kernel=best_kernel22)

predictions_GWR <- gwr.predict(
  formula = AKB ~   PK + RLS + KP+PKKP  + SL + AML +  IDL + ASI + BBLR + Rdokter + Rpuskesmas,
  data = spatial_train22,
  predictdata = spatial_test22,
  bw = best.bw22, adaptive = T,kernel=best_kernel22)

summary(best.gwr22$SDF$Local_R2)

# Hasil prediksi
FIPS22.test.xy$Pred.GWR <- predictions_GWR$SDF$prediction
FIPS22.train.xy$Pred.GWR <- best.gwr22$SDF$yhat
FIPS22.train.xy$Res.GWR <- spatial_train22$AKB - FIPS22.train.xy$Pred.GWR
FIPS22.test.xy$Res.GWR <- spatial_test22$AKB - FIPS22.test.xy$Pred.GWR
summary(FIPS22.test.xy$Res.GWR)

# Nilai target aktual dari data pelatihan
y_actual <- spatial_train22$AKB
y_pred <- FIPS22.train.xy$Pred.GWR
y_mean <- mean(y_actual)
TSS <- sum((y_actual - y_mean)^2)
RSS <- sum((y_actual - y_pred)^2)
R_squared_gwr<- 1 - (RSS / TSS)
R_squared_gwr

# Nilai target aktual dari data testing
y_actual_test <- spatial_test22$AKB
y_pred_test <- FIPS22.test.xy$Pred.GWR 
y_mean_test <- mean(y_actual_test)
TSS_test <- sum((y_actual_test - y_mean_test)^2)
RSS_test <- sum((y_actual_test - y_pred_test)^2)
R_squared_test <- 1 - (RSS_test / TSS_test)
R_squared_test


cat('GWR Test RMSE:', round(sqrt(mean((FIPS22.test.xy$Pred.GWR-FIPS22.test.xy$AKB)^2 , na.rm = TRUE)), digits=3), '\n')
cat('GWR Test MAE:', round(mean(abs(FIPS22.test.xy$Pred.GWR-FIPS22.test.xy$AKB) , na.rm = TRUE ), digits=3), '\n')
cat('GWR Test R2:', round(summary(lm(AKB~Pred.GWR,FIPS22.test.xy))$r.squared, digits=3), '\n')

cat('GWR RMSE:', round(sqrt(mean((FIPS22.train.xy$Pred.GWR-FIPS22.train.xy$AKB)^2 , na.rm = TRUE)), digits=3), '\n')
cat('GW MAE:', round(mean(abs(FIPS22.train.xy$Pred.GWR-FIPS22.train.xy$AKB) , na.rm = TRUE ), digits=3), '\n')
cat('GWR R2:', round(summary(lm(AKB~Pred.GWR,FIPS22.train.xy))$r.squared, digits=3), '\n')

#============= R2 Lokal GWR
#R2 local testing
y_actual_test <- FIPS22.test.xy$AKB
y_pred_test.gwr <- FIPS22.test.xy$Pred.GWR
y_mean_test.gwr <- mean(y_actual_test)
R_squared_local.gwr <- 1 - ((y_actual_test - y_pred_test.gwr)^2 / (y_actual_test - y_mean_test.gwr)^2)
R_squared_local.gwr

#R2 local training
y_actual_train <- FIPS22.train.xy$AKB
y_pred_train.gwr <- FIPS22.train.xy$Pred.GWR
y_mean_train.gwr <- mean(y_actual_train)
R_squared_local1.gwr <- 1 - ((y_actual_train - y_pred_train.gwr)^2 / (y_actual_train - y_mean_train.gwr)^2)
R_squared_local1.gwr

FIPS22.test.xy$R2.GWR <-  R_squared_local.gwr
FIPS22.train.xy$R2.GWR <-  R_squared_local1.gwr
FIPS22.test.xy$Standres.GWR <- ((FIPS22.test.xy$AKB-FIPS22.test.xy$Pred.GWR) - mean((FIPS22.test.xy$AKB-FIPS22.test.xy$Pred.GWR))) / sd((FIPS22.test.xy$AKB-FIPS22.test.xy$Pred.GWR))
FIPS22.train.xy$Standres.GWR <- (FIPS22.train.xy$Res.GWR - mean(FIPS22.train.xy$Res.GWR )) / sd(FIPS22.train.xy$Res.GWR)

# Membuat kategori untuk kolom R-Square Lokal
combined_data_GWRF <- combined_data_GWRF %>%
  mutate(
    R2.Local.GWR = case_when(
      R2.GWR <= 0.2 ~ "≤ 0.2",
      R2.GWR > 0.2 & R2.GWR <= 0.4 ~ "(0.2, 0.4]",
      R2.GWR > 0.4 & R2.GWR <= 0.6 ~ "(0.4, 0.6]",
      R2.GWR > 0.6 & R2.GWR <= 0.8 ~ "(0.6, 0.8]",
      R2.GWR > 0.8 ~ "> 0.8"
    )
  ) %>%
  # Convert the category column into a factor for efficiency
  mutate(R2.Local.GWR = factor(R2.Local.GWR, levels = c("≤ 0.2", "(0.2, 0.4]", "(0.4, 0.6]", "(0.6, 0.8]", "> 0.8")))

summary(combined_data_GWRF$R2.Local.GWR)

# Define colors using RColorBrewer palette
category_colors <- brewer.pal(n = 5, name = "Set3")  # Using Set3 palette with 5 colors

# Plot
r2.gwr22 <- ggplot(data = combined_data_GWRF) +
  geom_sf(aes(fill = R2.Local.GWR, color = dataset), show.legend = TRUE) +
  scale_fill_manual(values = category_colors, name = NULL) + # Warna berdasarkan kategori
  scale_color_manual(values = c("Testing" = "red"), name = NULL) +
  labs(title = NULL, x = "Longtitude", y = "Latitude") +
  theme_void() + # Menghilangkan grid dan background
  theme(
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
  ) +
  guides(
    fill = guide_legend(
      nrow = 1, byrow = TRUE,
      override.aes = list(color = category_colors) # Warna garis sesuai kategori
    ),
    color = guide_legend(
      override.aes = list(fill = "white", color = "red", size = 1.2) # Testing fill putih dengan garis merah
    )
  ) +
  coord_sf(xlim = c(105.5, 114.5), ylim = c(-8.8, -5.5)) # Zoom wilayah Jawa Timur

ggsave("D:/Kuliah S2/Thesis/Syntax/Eksplorasi/22_GWR_R2.png", plot = r2.gwr22, dpi = 800,width = 10, height = 5,  units = "in")

#============ MORAN LOKAL GWR FULL TRAINING TESTING
# Data Testing
coords22test <- test.df22[,c(4:5)]
neighbors22test <- knearneigh(coords22test, k = 5) 
weights22test <- nb2listw(knn2nb(neighbors22test), style = "W")

FIPS22.test.xy$Res.GWR <- FIPS22.test.xy$AKB - FIPS22.test.xy$Pred.GWR
FIPS22.test.xy$Standres.GWR <- (FIPS22.test.xy$Res.GWR - mean(FIPS22.test.xy$Res.GWR)) / sd(FIPS22.test.xy$Res.GWR)

# Moran's I dan Local Moran
moran_result.gwr <- moran.test(FIPS22.test.xy$Standres.GWR, weights22test)
print(moran_result.gwr)

moran_local_test.gwr <- localmoran(FIPS22.test.xy$Standres.GWR, weights22test)
FIPS22.test.xy$Local_Moran_I.gwr <- moran_local_test.gwr[, 1]
FIPS22.test.xy$p_value.gwr <- moran_local_test.gwr[, 5]

FIPS22.test.xy$Cluster.GWR <- with(FIPS22.test.xy, 
                                   ifelse(p_value.gwr < 0.05 & Res.GWR > 0 & Local_Moran_I.gwr > 0, "High-High",
                                          ifelse(p_value.gwr < 0.05 & Res.GWR > 0 & Local_Moran_I.gwr < 0, "High-Low",
                                                 ifelse(p_value.gwr < 0.05 & Res.GWR < 0 & Local_Moran_I.gwr > 0, "Low-Low",
                                                        ifelse(p_value.gwr < 0.05 & Res.GWR < 0 & Local_Moran_I.gwr < 0, "Low-High", "Not Significant")))))

FIPS22.test.xy$Cluster.GWR <- factor(FIPS22.test.xy$Cluster.GWR, 
                                     levels = c("High-High", "High-Low", "Low-Low", "Low-High", "Not Significant"))

factor_counts_test <- table(FIPS22.test.xy$Cluster.GWR)
factor_counts_test

# Data Training
FIPS22.train.xy$Res.GWR <- FIPS22.train.xy$AKB - FIPS22.train.xy$Pred.GWR
FIPS22.train.xy$Standres.GWR <- (FIPS22.train.xy$Res.GWR - mean(FIPS22.train.xy$Res.GWR)) / sd(FIPS22.train.xy$Res.GWR)

moran_result.gwrtrain <- moran.test(FIPS22.train.xy$Standres.GWR, weights22)
print(moran_result.gwrtrain)

moran_local.gwrtrain <- localmoran(FIPS22.train.xy$Standres.GWR, weights22)
FIPS22.train.xy$Local_Moran_I.gwr <- moran_local.gwrtrain [, 1]
FIPS22.train.xy$p_value.gwr <- moran_local.gwrtrain [, 5]
FIPS22.train.xy$Cluster.GWR <- with(FIPS22.train.xy, 
                                    ifelse(p_value.gwr < 0.05 & Res.GWR > 0 & Local_Moran_I.gwr > 0, "High-High",
                                           ifelse(p_value.gwr < 0.05 & Res.GWR > 0 & Local_Moran_I.gwr < 0, "High-Low",
                                                  ifelse(p_value.gwr < 0.05 & Res.GWR < 0 & Local_Moran_I.gwr > 0, "Low-Low",
                                                         ifelse(p_value.gwr < 0.05 & Res.GWR < 0 & Local_Moran_I.gwr < 0, "Low-High", "Not Significant")))))
FIPS22.train.xy$Cluster.GWR <- factor(FIPS22.train.xy$Cluster.GWR, 
                                      levels = c("High-High", "High-Low", "Low-Low", "Low-High", "Not Significant"))
factor_counts_train <- table(FIPS22.train.xy$Cluster.GWR)
factor_counts_train 

############## VISUALISASI LOCAL MORAN RESIDUAL FULL GWR
FIPS22.train.xy$dataset <- "Training"
FIPS22.test.xy$dataset <- "Testing"
colnames(FIPS22.test.xy)
colnames(FIPS22.train.xy)

# Gabungkan kedua data
combined_data_GWRF <- rbind(FIPS22.train.xy, FIPS22.test.xy)
combined_data_GWRF[combined_data_GWRF$Cluster.GWR == "High-High", ]
combined_data_GWRF[combined_data_GWRF$Cluster.GWR == "Low-Low", ]
combined_data_GWRF[combined_data_GWRF$Cluster.GWR == "High-Low", ]
combined_data_GWRF[combined_data_GWRF$Cluster.GWR == "Low-High", ]
nrow(combined_data_GWRF[combined_data_GWRF$Cluster.GWR == "Not Significant", ])
table(combined_data_GWRF$Cluster.GWR)

# Menentukan semua level untuk Cluster.GWR
combined_data_GWRF$Cluster.GWR <- factor(
  combined_data_GWRF$Cluster.GWR, 
  levels = c("High-High", "High-Low", "Low-Low", "Low-High", "Not Significant")
)

# Mendefinisikan warna untuk setiap cluster
cluster_colors <- c(
  "Not Significant" = "gray", 
  "Low-Low" = "lightblue", 
  "High-Low" = "pink", 
  "Low-High" = "blue", 
  "High-High" = "red"
)

# Membuat plot peta dengan ggplot
localmoran.gwr22 <- ggplot(data = combined_data_GWRF) +
  geom_sf(aes(fill = Cluster.GWR, color = dataset), show.legend = TRUE) +
  scale_fill_manual(values = cluster_colors, name = NULL, drop = FALSE) + # Warna berdasarkan kategori
  scale_color_manual(values = c("Testing" = "red"), name = NULL) +
  labs(title = "Local Moran I Residual GWR") +
  theme_void() + # Menghilangkan grid dan background
  theme(
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold")) +
  guides(
    fill = guide_legend(
      nrow = 1, byrow = TRUE,
      override.aes = list(color = "black")),
    color = guide_legend(
      override.aes = list(fill = "white", color = "red", size = 1.2) # Testing fill putih dengan garis merah
    )) +
  coord_sf(xlim = c(105.5, 114.5), ylim = c(-8.8, -5.5)) # Zoom wilayah Jawa Timur

# Simpan plot
ggsave("D:/Kuliah S2/Thesis/Syntax/Eksplorasi/22_GWR_local_moran_full.png", 
       plot = localmoran.gwr22, dpi = 800, width = 10, height = 5, units = "in")

#============ STANDAR RESIDUAL GWR
combined_data_GWRF$Category.GWR <- cut(combined_data_GWRF$Standres.GWR,
                                       breaks = 5,  # Membagi data menjadi 5 bagian
                                       #labels = c("≤ -1.05", "(-1.05,-0.106]", "(-0.106,0.834]", "(0.834,1.77]", "> 1.77"),
                                       include.lowest = TRUE)
combined_data_GWRF$Category.GWR <- as.factor(combined_data_GWRF$Category.GWR)
category_colorsres <- brewer.pal(n = 4, name = "Spectral")  # Using Set3 palette with 5 colors

# Membuat plot peta dengan ggplot
stanres.gwr <- ggplot(data = combined_data_GWRF) +
  geom_sf(aes(fill = Category.GWR, color = dataset), show.legend = TRUE) +
  scale_fill_manual(values = category_colorsres, name = NULL) +  # Menggunakan skala warna manual untuk kategori
  scale_color_manual(values = c("Testing" = "red"), name = NULL) +
  labs(title = "Residual Standarisasi") +
  theme_void() +  # Menghilangkan grid dan background
  theme(
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold")) +
  guides(
    fill = guide_legend(
      nrow = 1, byrow = TRUE,
      override.aes = list(color = category_colorsres)),
    color = guide_legend(
      override.aes = list(fill = "white", color = "red", size = 1.2))) +
  coord_sf(xlim = c(105.5, 114.5), ylim = c(-8.8, -5.5)) # Zoom wilayah Jawa Timur
ggsave("D:/Kuliah S2/Thesis/Syntax/Eksplorasi/22_GWR_res stand full.png", plot = stanres.gwr, dpi = 800, width = 10, height = 5, units = "in")

#====================== MODEL OLS
model.ols <- lm(AKB ~   PK + RLS + KP+PKKP  + SL + AML +  IDL + ASI + BBLR + Rdokter + Rpuskesmas  , data = spatial_train22)
summary(model.ols)

#===== PREDIKSI OLS
pred.ols <- predict(model.ols, spatial_test22)
pred.ols

FIPS22.test.xy$Pred.ols <- pred.ols
FIPS22.train.xy$Pred.ols <- predict(model.ols)
FIPS22.train.xy$Res.ols <- spatial_train22$AKB - FIPS22.train.xy$Pred.ols
FIPS22.test.xy$Res.ols <- spatial_test22$AKB - FIPS22.test.xy$Pred.ols

#============ EVALUASI MODEL
#==== DATA TRAINING
rmse.ols <- sqrt(mean((FIPS22.train.xy$Pred.ols-FIPS22.train.xy$AKB)^2))
mae.ols <- mean(abs(FIPS22.train.xy$Pred.ols-FIPS22.train.xy$AKB))

#==== DATA TESTING
cat('OLS Test RMSE:', round(sqrt(mean((FIPS22.test.xy$Pred.ols-FIPS22.test.xy$AKB)^2 , na.rm = TRUE)), digits=3), '\n')
cat('OLS Test MAE:', round(mean(abs(FIPS22.test.xy$Pred.ols-FIPS22.test.xy$AKB) , na.rm = TRUE ), digits=3), '\n')
cat('OLS Test R2:', round(summary(lm(AKB~Pred.ols,FIPS22.test.xy))$r.squared, digits=3), '\n')
