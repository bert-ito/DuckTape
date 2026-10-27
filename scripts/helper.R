defineFDR<-function(pred,xthreshold=0.05, upper=0.1, lower=0.01)
{
	totalp<-pred@tp[[1]]+pred@fp[[1]]
	fdr<-pred@fp[[1]]/totalp
	cutoffs<-pred@cutoffs[[1]]
	FDRinrange<- fdr[fdr<upper&fdr>lower&cutoffs<200] # this threshold needs further tuning.
	CUTOFFinrange<-cutoffs[fdr<upper&fdr>lower&cutoffs<200]
	lo<-loess(CUTOFFinrange ~ FDRinrange)
	thr<-predict(lo,newdata=c(xthreshold))
	return(thr)
}


gtf2genePos <- function(myfile="test.gtf") {
	# this is a helper function to get gene positions from GTF files
	# works for Arabidopsis GTF file,
	# for other organisms, change the chronly variable.
	gtf <- read.delim(myfile, header=FALSE)
	colnames(gtf) <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame",     
			"attributes")
	chronly <- c("Chr1","Chr2","Chr3","Chr4","Chr5","ChrM","ChrC")
	gtf <- gtf[as.character(gtf$seqname) %in% chronly, ] # Cleanup to remove non-chromosome rows
	gtf[, length(gtf[1,])] <- gsub(" {1,}|gene_id|;.*", "", gtf$attributes)
	gtf1 <- gtf[order(gtf$seqname, gtf$start), ]
	gtf1 <- gtf1[!duplicated(gtf1$attributes), ] # Returns only rows with most left gene positions.
	gtf2 <- gtf[order(gtf$seqname, -gtf$end), ]
	gtf2 <- gtf2[!duplicated(gtf2$attributes), ] # Returns only rows with most right gene positions.
	gtf1 <- gtf1[order(gtf1$attributes), ]
	gtf2 <- gtf2[order(gtf2$attributes), ]
	gtf1[,5] <- gtf2[,5] # Assigns proper end postions
	gtfgenepos <- gtf1 # Returns gtf containing only start/end positions for all genes
	gtfgenepos <- gtfgenepos[order(gtfgenepos$seqname, gtfgenepos$start), ]
	return(gtfgenepos)
}

getnegativepairs<-function(gtf, dist = 500)
{
	ngenes<-nrow(gtf)
	leftgene<- gtf[1:(ngenes-3),c('seqname','start','end','strand','attributes')]
	midgene1<- gtf[2:(ngenes-2),c('seqname','start','end','strand','attributes')]
	midgene2<- gtf[3:(ngenes-1),c('seqname','start','end','strand','attributes')]
	rightgene<- gtf[4:ngenes,c('seqname','start','end','strand','attributes')]
	bothgenes<-cbind(leftgene, midgene1, midgene2, rightgene)
	midleftgenestrand<-midgene1[,'strand']
	midrightgenestrand<-midgene2[,'strand']
	out<-bothgenes[(midleftgenestrand=="+"&midrightgenestrand=="-")|
					(midleftgenestrand=="-"&midrightgenestrand=="+"),]
	out<-out[out[,6]==out[,11],]
	out<-out[(out[,7]-out[,3])>=dist&(out[,12]-out[,8])>=dist&(out[,17]-out[,13])>=dist,c(10,15,3,7,8,12,13,17)]
	rownames(out)<-NULL
	colnames(out)<-c('locus1','locus2')
	return(out[,c(1,2)])
}



getBICmatsol<-function(smat,asmat,labels,llanti_llanti_full,llsense=llsense_full)
{
	# this version get not only the BIC, but also the hill climbing solution of the read count
	stopifnot(ncol(smat)==length(labels),ncol(asmat)==length(labels), nrow(smat)==nrow(asmat))
	allmat<-smat+asmat;	conds<-as.factor(labels)
	nreplicates<-ncol(asmat);
	cds<-newCountDataSet(allmat,conds);	cds<-estimateSizeFactors(cds);	cds<-estimateDispersions(cds) 
	SF<-sizeFactors( cds );	rawVarFunc<-fitInfo(cds)$dispFun
	smatsum<-apply(smat,2,sum); asmatsum<-apply(asmat,2,sum); p<-asmatsum/(smatsum+asmatsum)
	
	ngenes <- nrow(smat); nconds <- length(labels)
	BICcormat<-rep(0.0,ngenes*(8+2*nreplicates))
	dim(BICcormat)<-c(ngenes,8+2*nreplicates)
	
	#predns: predicted number of reads from the sense gene
	#prednas: predicted number of reads from the antisense gene
	
	# this is not a flexable way of coding the column names.
	colnames(BICcormat)<-c("BIC","AveN1","AveN2",
			"p","preda1","predb1",
			"predns","prednas",
			'predns1','predns2','predns3',
			'prednas1','prednas2','prednas3')	
	
	rownames(BICcormat)<-rownames(smat)
	for(i in 1:ngenes)
	{
		N1<-smat[i,]
		N2<-asmat[i,]
		tmpx1<-c(floor(N1*p),floor(N2*p))
		tmpx2<-c(N2,N1)
		startx<-apply(rbind(tmpx1,tmpx2),2,min)
		# filter out small counts
		if(sum((N1+N2)<5)>1)
		{
			BICcormat[i,1]<-c(-Inf)
			next
		}
		
		sol<-hillclimb(startx,rep(0,nconds*2),c(N2,N1),llanti,N1,N2,p,SF,rawVarFunc)
		out<-BICanti(sol,c(N1,N2),p,SF,rawVarFunc,llanti,llsense)
		BICcormat[i,1]<-round(out,2)
		BICcormat[i,2]<-round(mean(N1),2)
		BICcormat[i,3]<-round(mean(N2),2)
		BICcormat[i,4]<-round(mean(N2/(N1+N2)),2)
		
		nsample = length(labels)
		BICcormat[i,5]<-round(mean(sol[1:nsample]),2)
		BICcormat[i,6]<-round(mean(sol[(nsample+1):(nsample+3)]),2)
		
		predns <- N1 - sol[(nsample+1):(nsample+nsample)] + sol[1:nsample]
		prednas <- N2 - sol[1:nsample] + sol[(nsample+1):(nsample+nsample)]
		BICcormat[i,7]<-round(mean(predns),2)
		BICcormat[i,8]<-round(mean(prednas),2)
		BICcormat[i,9:11]<-predns
		BICcormat[i,12:14]<-prednas
	}
	names(BICcormat)<-rownames(smat)[1:ngenes]
	BICor <- order(BICcormat[,1],decreasing=T)
	BICcormat<-BICcormat[BICor,]
	return(BICcormat)
}

getNASTIscore<-function(smat,asmat, llanti=llanti_full, llsense=llsense_full)
{
	labels<-rep('wr',ncol(smat))
	stopifnot(ncol(smat)==length(labels),ncol(asmat)==length(labels), nrow(smat)==nrow(asmat))
	allmat<-smat+asmat; conds<-as.factor(labels)
	cds<-newCountDataSet(allmat,conds); cds<-estimateSizeFactors(cds); cds<-estimateDispersions(cds) 
	SF<-sizeFactors( cds );	rawVarFunc<-fitInfo(cds)$dispFun
	smatsum<-apply(smat,2,sum); asmatsum<-apply(asmat,2,sum)
	p<-asmatsum/(smatsum+asmatsum)
	
	ngenes <- nrow(smat); nconds <- length(labels)
	BICcormat<-rep(0.0,ngenes)
	for(i in 1:ngenes)
	{
		N1<-smat[i,]
		N2<-asmat[i,]
		tmpx1<-c(floor(N1*p),floor(N2*p))
		tmpx2<-c(N2,N1)
		startx<-apply(rbind(tmpx1,tmpx2),2,min)
		if(sum((N1+N2)<5)>1){
			BICcormat[i]<-c(-Inf)
			next}
		
		sol<-hillclimb(startx,rep(0,nconds*2),c(N2,N1),llanti,N1,N2,p,SF,rawVarFunc)
		out<-BICanti(sol,c(N1,N2),p,SF,rawVarFunc,llanti,llsense)
		BICcormat[i]<-out
	}
	names(BICcormat)<-rownames(smat)[1:ngenes]
	BICcormat<-sort(BICcormat,decreasing=T)
	return(BICcormat)
}

NASTIpredict<-function(poscounts, negcounts, pospairs, negpairs, llanti=llanti_full, llsense=llsense_full)
{
	allBIC<-getNASTIscore(poscounts,negcounts)
	bicout<-EvalBIC(allBIC,poscounts,negcounts,pospairs,negpairs)
	out<-list()
	out$predictions<-bicout[,1]
	out$labels<-bicout[,2]
	return(out)
}

EvalBIC<-function(bicmat, smat, asmat, WangPairs, noneighpairs, mincount=5)
{
	
	allcounts<-smat+asmat
	rownames(allcounts)<-rownames(smat)
	nconds <- ncol(allcounts)
	
	WangLeft<-WangPairs[,1]; WangRight<-WangPairs[,2]
	WangLeftExp<-allcounts[WangLeft,]; WangRightExp<-allcounts[WangRight,]
	WangLeftExp[WangLeftExp<mincount]=0; WangLeftExp[WangLeftExp>=mincount]=1
	WangRightExp[WangRightExp<mincount]=0; WangRightExp[WangRightExp>=mincount]=1
	
	tmpLeft<-apply(WangLeftExp[,c(1:nconds)],1,sum); tmpRight<-apply(WangRightExp[,c(1:nconds)],1,sum)
	WangLeftKeep<-names(tmpLeft[tmpLeft==nconds]);	WangRightKeep<-names(tmpRight[tmpRight==nconds])
	WangKeepPairs<- WangPairs[WangPairs[,1]%in%WangLeftKeep&WangPairs[,2]%in%WangRightKeep,]
	
	noneighleft<-noneighpairs[,1]; noneighright<-noneighpairs[,2]
	noneighleftExp<-allcounts[noneighleft,]; noneighrightExp<-allcounts[noneighright,]
	
	noneighleftExp[noneighleftExp<mincount]=0; noneighleftExp[noneighleftExp>=mincount]=1
	noneighrightExp[noneighrightExp<mincount]=0; noneighrightExp[noneighrightExp>=mincount]=1
	
	tmpLeft<-apply(noneighleftExp[,c(1:nconds)],1,sum); tmpRight<-apply(noneighrightExp[,c(1:nconds)],1,sum)
	noneighleftkeep<-names(tmpLeft[tmpLeft==nconds]); noneighrightkeep<-names(tmpRight[tmpRight==nconds])
	noneighKeepPairs<- noneighpairs[noneighpairs[,1]%in%noneighleftkeep&noneighpairs[,2]%in%noneighrightkeep,c(1,2)]
	
	PosBIC<-cbind(bicmat[WangKeepPairs[,1]],bicmat[WangKeepPairs[,2]])
	outBICpos<-apply(PosBIC,1,max)
	nposbic <- length(outBICpos)
	
	NegBIC<-cbind(bicmat[noneighKeepPairs[,1]],bicmat[noneighKeepPairs[,2]])
	outBICneg<-apply(NegBIC,1,max)
	outBICneg<-bicmat[noneighKeepPairs[,1]]
	outBIC<-c(outBICpos,outBICneg)
	outlabel<-c(rep(1,length(outBICpos)),rep(0,length(outBICneg)))
	return(cbind(outBIC,outlabel))
}

FindNATs<-function(allBIC, tBIC, filterpairs, annotation=NULL, distance = 500)
{
	
	outgenes<-allBIC[allBIC>tBIC]
	
	#remove genes that are in the filterpairs 
	tmp<-outgenes[!names(outgenes)%in%filterpairs[,1]]
	tmp1<-tmp[!names(tmp)%in%filterpairs[,2]] 
	
	#known_pos_genes are genes that are from the positive training sets.
	known_pos_left_genes<-names(outgenes[names(outgenes)%in%filterpairs[,1]])
	known_pos_right_genes<-names(outgenes[names(outgenes)%in%filterpairs[,2]])
	known_pos_genes<-unique(c(known_pos_left_genes,known_pos_right_genes))
	
	#keep genes that passed selection
	firsttokeep<-filterpairs[filterpairs[,1]%in%names(outgenes),]
	sectokeep<-filterpairs[filterpairs[,2]%in%names(outgenes),]
	
	known_pos_genepairs<-rbind(firsttokeep,sectokeep)
	
	#newasgenes are those genes that passed the threshold,
	#next is to match the neighbor gene that potentially has the correct orientation
	newasgenes<-sort(tmp1,decreasing=T)
	
	rownames(annotation)<-annotation[,9]
	tmp<-annotation[!annotation[,2]%in%c('snRNA','tRNA','miRNA','snoRNA'),] # filter out small RNAs
	annotation<-tmp
	
	# keep those appear in annotation file
	tmp<-newasgenes[names(newasgenes)%in%rownames(annotation)]
	newasgenes<-tmp
	asgeneannotation<-annotation[names(newasgenes),]
	
	Ngenetocheck = 5
	
	allcandidates<-list()
	for (eachgene in names(newasgenes))
	{
		geneloc <- which(rownames(annotation)==eachgene)

		# if gene is located at the end of annotation files.
		geneleftloc <- geneloc-Ngenetocheck
		if (geneleftloc<0) geneleftloc<-1
		generightloc <- geneloc+Ngenetocheck
		if (generightloc>nrow(annotation)) generightloc<-nrow(annotation)
		
		surroundinggenes<-annotation[geneleftloc:generightloc,]
		
		#for the current gene, get the strand,left and right position
		eachgenestrand<-surroundinggenes[eachgene,7]
		eachgeneleft<-surroundinggenes[eachgene,4]
		eachgenerigth<-surroundinggenes[eachgene,5]
		
		#check each gene in the list of neighbor genes	
		candidates<-NULL
		for (eachnextgene in rownames(surroundinggenes))
		{
			if (eachnextgene == eachgene) next
			if (surroundinggenes[eachnextgene,7] == eachgenestrand) next
			
			leftdist<-c(eachgeneleft-surroundinggenes[eachnextgene,5])
			rightdist<-c(surroundinggenes[eachnextgene,4]-eachgenerigth)
			
			if (leftdist>distance | rightdist>distance) next
			candidates<-c(candidates,eachnextgene)
		}
		
		allcandidates[[ eachgene ]]<-candidates
	}
	
	#outpairs hold the best pairs for each predicted antisense gene
	#if a gene has no pair, then return NA
	
	outpairs<-data.frame(cbind(rep('',nrow(asgeneannotation)),
					rep('',nrow(asgeneannotation))),
			stringsAsFactors=FALSE)
	
	rownames(outpairs)<-rownames(asgeneannotation)
	colnames(outpairs)<-c('genename','mark')
	outgenenames<-rownames(asgeneannotation)	
	for (eachgene in outgenenames)
	{
		if (is.null(allcandidates[[eachgene]]))
		{
			#no gene the target gene
			outpairs[eachgene,1]="NOANNOT"
			outpairs[eachgene,2]=paste(">",distance,sep='')
		}
		else if(length(allcandidates[[eachgene]])==1)
		{
			#For genes that only has one candidate neighbor within a given distance.
			neighborgene<-allcandidates[[eachgene]]
			neighborannot<-annotation[neighborgene,]
			x1<-asgeneannotation[eachgene,4]-neighborannot[,4]
			x2<-asgeneannotation[eachgene,4]-neighborannot[,5]
			x3<-asgeneannotation[eachgene,5]-neighborannot[,4]
			x4<-asgeneannotation[eachgene,5]-neighborannot[,5]
			if ((sign(x1)>0&sign(x2)<=0&sign(x3)>0&sign(x4)>=0)|
					(sign(x1)<=0&sign(x2)<0&sign(x3)>0&sign(x4)>=0)|
					(sign(x1)<=0&sign(x2)<0&sign(x3)>=0&sign(x4)<0)|
					(sign(x1)>0&sign(x2)<0&sign(x3)>0&sign(x4)<0))
			{
				#four possible antisense patterns:
				#+(-0)+(+0), (-0)-+(+0), (-0)-(+0)-, +-+-
				outpairs[eachgene,1]=neighborgene
				outpairs[eachgene,2]="ANNOTATION"
			}
			else
			{
				outpairs[eachgene,1]=neighborgene
				outpairs[eachgene,2]="NOTANNOT"
			}
		}
		else
		{
			neighborgene<-allcandidates[[eachgene]]
			neighborannot<-annotation[neighborgene,]
			x1<-asgeneannotation[eachgene,4]-neighborannot[,4]
			x2<-asgeneannotation[eachgene,4]-neighborannot[,5]
			x3<-asgeneannotation[eachgene,5]-neighborannot[,4]
			x4<-asgeneannotation[eachgene,5]-neighborannot[,5]
			
			#check each candidate for genes that are annotated as overlapping genes.
			# this step will find the first annotated gene that could explain the antisense transcription
			# however, there could be cases where more than one gene are possible candidate.
			foundannotation<-FALSE
			for (eachi in 1:length(neighborgene)){
				if ((sign(x1[eachi])>0&sign(x2[eachi])<=0&sign(x3[eachi])>0&sign(x4[eachi])>=0)|
						(sign(x1[eachi])<=0&sign(x2[eachi])<0&sign(x3[eachi])>0&sign(x4[eachi])>=0)|
						(sign(x1[eachi])<=0&sign(x2[eachi])<0&sign(x3[eachi])>=0&sign(x4[eachi])<0)|
						(sign(x1[eachi])>0&sign(x2[eachi])<0&sign(x3[eachi])>0&sign(x4[eachi])<0))
				{
					outpairs[eachgene,1]=neighborgene[eachi]
					outpairs[eachgene,2]="ANNOTATION"
					foundannotation<-TRUE
					break
				}
				else
				{
					next
				}
			}
			if (foundannotation) next
			
			
			x1<-abs(asgeneannotation[eachgene,4]-neighborannot[,4])
			x2<-abs(asgeneannotation[eachgene,4]-neighborannot[,5])
			x3<-abs(asgeneannotation[eachgene,5]-neighborannot[,4])
			x4<-abs(asgeneannotation[eachgene,5]-neighborannot[,5])
			
			#center difference
			x5<-abs((asgeneannotation[eachgene,5]+asgeneannotation[eachgene,4])/2
							-(neighborannot[,5]+neighborannot[,4])/2)
			
			disttoasgene<-apply(cbind(x1,x2,x3,x4,x5),1,min)
			names(disttoasgene)<-neighborgene
			#print(cbind(x1,x2,x3,x4,x5))
			if(eachgene == 'AT2G36360'){
				print(disttoasgene)
				print(neighborannot)
				print(asgeneannotation[eachgene,])
				print(cbind(x1,x2,x3,x4,x5))
				return(NULL)
			}
			disttoasgene<-sort(disttoasgene)
			
			closestgene <- names(disttoasgene)[1]
			
			neighborgene<-closestgene
			neighborannot<-annotation[neighborgene,]
			x1<-asgeneannotation[eachgene,4]-neighborannot[,4]
			x2<-asgeneannotation[eachgene,4]-neighborannot[,5]
			x3<-asgeneannotation[eachgene,5]-neighborannot[,4]
			x4<-asgeneannotation[eachgene,5]-neighborannot[,5]
			if ((sign(x1)>0&sign(x2)<=0&sign(x3)>0&sign(x4)>=0)|
					(sign(x1)<=0&sign(x2)<0&sign(x3)>0&sign(x4)>=0)|
					(sign(x1)<=0&sign(x2)<0&sign(x3)>=0&sign(x4)<0)|
					(sign(x1)>0&sign(x2)<0&sign(x3)>0&sign(x4)<0))
			{
				#four possible antisense patterns:
				#+(-0)+(+0), (-0)-+(+0), (-0)-(+0)-, +-+-
				outpairs[eachgene,1]=neighborgene
				outpairs[eachgene,2]="ANNOTATION"
			}
			else
			{
				outpairs[eachgene,1]=neighborgene
				outpairs[eachgene,2]="NOTANNOT"
			}
		}
	}
	
	#combine genes with known pairs
	byTAIR<-outpairs[outpairs[,2]=="ANNOTATION",]
	byTAIR[,2]<-rownames(byTAIR)
	colnames(known_pos_genepairs)<-colnames(byTAIR)
	known_pos_genepairs<-rbind(known_pos_genepairs,byTAIR)
	
	#for known pairs, filter those that appeared twice 
	outknownpairstmp <- NULL
	for (i in 1:nrow(known_pos_genepairs))
	{
		tmpname<-paste(sort(c(known_pos_genepairs[i,1],known_pos_genepairs[i,2])),collapse='_')
		outknownpairstmp<- c(outknownpairstmp,tmpname)
	}
	outknownpairstmp <- unique(outknownpairstmp)
	for (i in 1:length(outknownpairstmp))
	{
		if (i == 1)
		{
			outknownpairsfilter<-strsplit(outknownpairstmp[i],'_')[[1]]
		}
		else
		{
			outknownpairsfilter<-rbind(outknownpairsfilter,strsplit(outknownpairstmp[i],'_')[[1]])
		}
	}
	
	newpairs<-outpairs[outpairs[,2]=='NOTANNOT',]
	
	#filter out genes that are potentially be antisense to two different genes.
	#there are four such pairs in WT prediction.
	# so that the number is small and not a big problem for subsequent analysis.
	outnewpairs<-cbind(newpairs[,1],rownames(newpairs))
	dupgenes<-names(which(table(newpairs[,1])>1))
	outnewpairs<-outnewpairs[!outnewpairs[,1]%in%dupgenes,]
	
	#also keep only one pair.
	outpairstmp <- NULL
	for (i in 1:nrow(outnewpairs))
	{
		tmpname<-paste(sort(c(outnewpairs[i,1],outnewpairs[i,2])),collapse='_')
		outpairstmp <- c(outpairstmp,tmpname)
	}
	outpairstmp<-unique(outpairstmp)
	outnewpairsfilter<-NULL
	for (i in 1:length(outpairstmp))
	{
		if (i == 1)
		{
			outnewpairsfilter<-strsplit(outpairstmp[i],'_')[[1]]
		}
		else
		{
			outnewpairsfilter<-rbind(outnewpairsfilter,strsplit(outpairstmp[i],'_')[[1]])
		}
	}
	
	
	neworphan<-outpairs[outpairs[,2]=='>500',]
	
	result<-NULL
	result$knownpairs<-outknownpairsfilter # training data and genes from TAIR.
	rownames(result$knownpairs)<-NULL
	result$newpairs<-outnewpairsfilter #newly found gene pairs
	rownames(result$newpairs)<-NULL
	result$neworphan<-rownames(neworphan) #orphan genes

	return(result)
	
	
}


