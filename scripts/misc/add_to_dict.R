dict$age <- rep(NA,nrow(dict))
dict$os <- rep(NA,nrow(dict))
dict$ptv80 <- rep(NA,nrow(dict))


for(a in levels(dict$patient)){
	dict$age[which(dict$patient==a)] <- tmp$age_at_diagnosis[which(tmp$patient==a)] 
	dict$os[which(dict$patient==a)] <- tmp$os[which(tmp$patient==a)] 
	dict$ptv80[which(dict$patient==a)] <- tmp$PTV80[which(tmp$patient==a)] 
}
