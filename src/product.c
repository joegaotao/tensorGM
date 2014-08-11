//byDim index from 1 on which product
void product(double *a, double *b, int *dimA, int *dimB, int *dimT,double *c,int *byDim)
{
	int sizeA;
	int i,j;
	int by=byDim[0]-1;
	int dimTotal=dimT[0];
	
	sizeA=1;
	for(i=0;i<dimTotal;i++) sizeA=sizeA*dimA[i];
	
	
	int* mappingC; // mapping coordinate to pos for matrix C
	mappingC=(int *) malloc(dimTotal*sizeof(int));
	mappingC[dimTotal-1]=1;
	for(i=dimTotal-2;i>=0;i--){
		if(i!=by-1){
			mappingC[i]=mappingC[i+1]*dimA[i+1];
		}else{
			mappingC[i]=mappingC[i+1]*dimB[0];
		}
	}
	
	int* coord;
	int pos;
	coord=(int *) malloc(dimTotal*sizeof(int));
	for(i=0;i<dimTotal;i++){
		coord[i]=0;
	}
	coord[dimTotal-1]=-1;
	//main computation
	for(i=0;i<sizeA;i++){
		coord[dimTotal-1]=coord[dimTotal-1]+1;
		for(j=dimTotal-1;j>0;j--){
			if(coord[j]==dimA[j]){
				coord[j]=0;
				coord[j-1]++;
			}
		}
		pos=0;
		for(j=0;j<dimTotal;j++){
			if(j!=by) pos+=coord[j]*mappingC[j];
		}
		for(j=0;j<dimB[0];j++){
			c[pos]+=a[i]*b[j*dimB[1]+coord[by]];
			pos+=mappingC[by];
		}
	}
}