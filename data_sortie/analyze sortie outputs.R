files = dir('outputs',pattern = '*out',full.names = TRUE)

i=1
df_this = read.table(files[i], sep = '\t', skip=5, header=TRUE)

#df_this[nrow(df_this),grep("Sdl.Abs.Den.",names(df_this),fixed=TRUE)]
n_final_this = tail(df_this[,grep("Sapl.Abs.Den.",names(df_this),fixed=TRUE)],1) + tail(df_this[,grep("Adult.Abs.Den.",names(df_this),fixed=TRUE)],1)

n_last_50yr_this = tail(df_this[,grep("Sapl.Abs.Den.",names(df_this),fixed=TRUE)],10) + tail(df_this[,grep("Adult.Abs.Den.",names(df_this),fixed=TRUE)],10)

n_cv_last_50yr_this = apply(n_last_50yr_this, 2, function(x) { cv = sd(x)/mean(x); cv[is.nan(cv)]=0; cv   })

n_cv_last_50yr_mean_this = mean(n_cv_last_50yr_this)

fn_this = gsub("_[0-9]\\.out","",gsub("GMF_","",basename(files[i])))
