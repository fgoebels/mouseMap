
library(GEOquery)
library(SRAdb)

geoIDs <- c(Meta(getGEO("GPL9250"))$sample_id, Meta(getGEO("GPL11002"))$sample_id, Meta(getGEO("GPL13112"))$sample_id, Meta(getGEO("GPL17021"))$sample_id)

valGEOIDs = unlist(lapply(geoIDs, getrnaseqGEOIDs))
output = unlist(lapply(valGEOIDs, printGEO))

getrnaseqGEOIDs <- function(x, geodir="~/workspace/GEO/GSM/"){
	
	geodata=Meta(getGEO(x, destdir=geodir))
	if(isNULL(geodata$library_source)){return(NULL)}
	if(isNULL(geodata$library_selection)){return(NULL)}
	if(isNULL(geodata$library_strategy)){return(NULL)}
	if(isNULL(geodata$organism_ch1)){return(NULL)}
	if(geodata$library_source=="transcriptomic" & geodata$library_selection=="cDNA" & geodata$library_strategy=="RNA-Seq" & geodata$organism_ch1=="Mus musculus"){
		return(x)
	}
}

isNULL <- function(x){
	return(typeof(x)=="NULL")
}

printGEO <- function(x, geodir="~/workspace/GEO/GSM/"){
	geodata=Meta(getGEO(x, destdir=geodir))
	rels = geodata$relation
	srxIDs = unlist(lapply(rels, getSRXID))
	srrIDs = unlist(lapply(srxIDs, mapSRXID))
	if(length(srrIDs)==0 | length(srrIDs)>1){
		return(NULL)
	}
	srxIDs = strjoin(srxIDs)
	srrIDs = strjoin(srrIDs)
	character = strjoin(geodata$characteristics_ch1)
	series_id = strjoin(geodata$series_id)
	out = strjoin(c(x, series_id, geodata$source_name_ch1, geodata$platform_id, geodata$instrument_model, geodata$molecule_ch1, character, srxIDs, srrIDs), sep="\t")
	return(out)
}

strjoin <- function(x, sep = ";"){
	return(paste(x, collapse=sep))
}

dwldGEOfiles <- function(x, geodir="~/workspace/GEO/GSM/"){
	out = NULL
	if(!file.exists(paste(c(geodir, x, ".soft"), collapse=""))){
		out = tryCatch(getGEOfile(x, destdir=geodir, amount="brief"), error=function(e) FALSE)
	}
	return(out) 
}

getSRXID <- function(x){
	if(substr(x,0,4)=="SRA:"){
		return(strsplit(x, "=")[[1]][-1])
	}
}


mapSRXID <- function(x, sqlFile = "/Users/florian/SRAmetadb.sqlite"){
	sra_con <- dbConnect(SQLite(), sqlFile)
	return(sraConvert(in_acc=x, sra_con=sra_con)$run)
}
