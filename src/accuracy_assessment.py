@author David T. Milodowski
@email d.t.milodowski@ed.ac.uk
@date 28/04/2016

import numpy as np

"""
This code contains the accuracy assessment routines outlined in the Olofsson et al., Remote Sensing of Environment 148 (2014) 42-57.

There is code to get sample sizes for a stratified random sample procedure.  There is also code to calculate the accuracy statistics of a confusion matrix, in addition to code to calculate the accuracy statistics for a stratified random sample accuracy assessment.
"""

"""
The first set of functions derive appropriate sample sizes, and perform the stratified random sample
"""

def get_sample_sizes(MappedClassAreas,PredictedClassUA,TargetSErr=0.01, n_min=75.):
    W=MappedClassAreas/np.sum(MappedClassAreas)
    S=np.sqrt(PredictedClassUA*(1.-PredictedClassUA))
    n = np.round(np.sum(W*S/TargetSErr)**2)
    print W
    ClassSampleSizes_iter = np.round(n*W)
    n_remain = n
#    if np.min(ClassSampleSizes_init)<n_min:
    while np.sum(ClassSampleSizes_iter<n_min)!=0:
        #count number of classes with smaller sample sizes
        N_small = np.sum(ClassSampleSizes_iter<=n_min)
        N_large = np.sum(ClassSampleSizes_iter>n_min)
        n_remain = n_remain-N_small*n_min
        
        W[ClassSampleSizes_iter<=n_min]=0
        W=W/np.sum(W)
        ClassSampleSizes_iter=np.round(n_remain*W)
        ClassSampleSizes_iter[ClassSampleSizes_iter<=n_min]=n_min

        if n_remain<0:
            print "NEED BIGGER SAMPLE SIZE"
            ClassSampleSizes_iter[:]=n_min

    ClassSampleSizes=ClassSampleSizes_iter.copy()
    print ClassSampleSizes
    return n, ClassSampleSizes

def retrieve_stratified_random_sample(ClassMap,ClassKeys,PredictedClassUA,XMin,YMax,XResolution,YResolution,TargetSErr=0.01, n_min=75.):
    
    N_classes = ClassKeys.size
    # create row and column id's
    nrows,ncols = ClassMap.shape
    rows_matrix = np.zeros((nrows,ncols))
    cols_matrix = np.zeros((nrows,ncols))
    for i in range(0,nrows):
        rows_matrix[i,:]=i
    for j in range(0,ncols):
        cols_matrix[:,j]=j

    MappedClassSizes = np.zeros(N_classes)
    for kk in range(0,N_classes):
        MappedClassSizes = np.sum(ClassMap==ClassKeys[kk])

    temp, ClassSampleSizes = get_sample_sizes(MappedClassSizes,PredictedClassUA,TargetSErr,n_min)

    N_samples = np.sum(ClassSampleSizes)
    sample_points = np.zeros((N_samples,5)) # this will be an array of points with values x,y,row,col,class

    count = 0
    
    for kk in range(0,N_classes):
        class_rows = rows_matrix[ClassMap==ClassKeys[kk]]
        class_cols = cols_matrix[ClassMap==ClassKeys[kk]]
        indices = np.arange(0,MappedClassAreas[kk])
        
        #randomly sample indices without replacement
        sample_indices = np.random.choice(indices, ClassSampleSizes[kk], replace=False)
        
        for ii in range(0,ClassSampleSizes[kk]):
            sample_points[ii,0] = class_cols[sample_indices[ii]]*XResolution + Xmin + XResolution/2 #x
            sample_points[ii,1] = Ymax - YResolution/2 - class_rows[sample_indices[ii]]*Y_Resolution #y
            sample_points[ii,2] = class_rows[sample_indices[ii]]#row
            sample_points[ii,3] = class_cols[sample_indices[ii]]#col
            sample_points[ii,4] = ClassKeys[kk]#class

    return sample_points

def write_sample_points_to_csv(sample_points,FILENAME,SAVEDIR='./'):
    hdr = "x,y,row,col,class"
    np.savetxt(SAVEDIR+FILENAME+'.csv',sample_points,",",header=hdr)


"""
The following functions calculate the accuracy statistics and area changes, including confidence intervals as outlined by Olofsson et al., 2014.
"""
def calculate_accuracy_stats(ConfusionMatrix):

    N_points = np.sum(ConfusionMatrix)
    q,temp = ConfusionMatrix.shape
    OA = np.sum(ConfusionMatrix.diagonal())/np.sum(ConfusionMatrix)
    UA = np.zeros(q)
    for i in range(0,q):
        UA[i]=ConfusionMatrix[i,i]/np.sum(ConfusionMatrix[i,:])

    PA = np.zeros(q)
    for i in range(0,q):
        PA[i]=ConfusionMatrix[i,i]/np.sum(ConfusionMatrix[:,i])
    
    return OA, UA, PA

def calculate_accuracy_stats_from_sample(CM, MappedClassAreas):

    q,temp = CM.shape
    n=np.sum(CM,axis=1)
 
    # First convert confusion matrix from representing pixels to representing proportional area
    CM_Estimator = CM.copy()
    W=MappedClassAreas/np.sum(MappedClassAreas)
    for i in range(0,q):
        for j in range(0,q):
            CM_Estimator[i,j]=W[i]*CM[i,j]/n[i]

    OA,UA,PA = calculate_accuracy_stats(CM_Estimator)

    # Now get sampling variance
    # i) Overall Accuracy
    Var_OA = np.sum(W*W*UA*(1-UA)/(n-1))
    # ii) User's accuracy
    Var_UA = UA*(1-UA)/(n-1)
    # iii) Producer's accuracy
    N_i=np.sum(CM_Estimator,axis=1)
    N_j=np.zeros(q)
    for j in range(0,q):
        for i in range(0,q):
            N_j[j]+=(N_i[i]/n[i])*CM[i,j]
    Term2 = np.zeros(q)
    for j in range(0,q):
        for i in range(0,q):
            if i!=j:
                Term2[j] += N_i[i]**2 * CM[i,j]/n[i] * (1-CM[i,j]/n[i]) / (n[i]-1)

    Var_PA =(1/N_j**2) * ( (N_i**2) * ((1-PA)**2) * UA * (1-UA) / (n-1) + (PA**2) * (Term2) )
    
    # Now get Area estimate for stratified estimator (as proportion of whole area)
    p_area=np.zeros(q)
    for k in range(0,q):
        for i in range(0,q):
            p_area[k]+=W[i]*CM[i,k]/n[i]

    # Now get standard error for this estimate (as proportion)
    Serr = np.zeros(q)
    for k in range(0,q):
        for i in range(0,q):
            Serr[k]+=(W[i]*CM_Estimator[i,k]-CM_Estimator[i,k]**2) / (n[i]-1)
    Serr_area=np.sqrt(Serr)
 
    """
    print n
    print "N_j"
    print N_j

    print "term 1"
    print  (N_i**2) * ((1-PA)**2) * UA * (1-UA) / (n-1)

    print "term 2"
    print (PA**2) * (Term2)

    print "V(P)"
    print Var_PA

    print "Confidence Ints"
    print "OA"
    print 1.96*np.sqrt(Var_OA)
    print "UA"
    print 1.96*np.sqrt(Var_UA)
    print "PA"
    print 1.96*np.sqrt(Var_PA)
    print " "

    print p_area
    print Serr_area
    """
    return OA, UA, PA, Var_OA, Var_UA, Var_PA, p_area, Serr_area

#MappedClassAreas = np.array([200000., 150000., 3200000., 6450000.])
