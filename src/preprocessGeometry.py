# Written by Markus Muhr @ TUM 2021

import numpy as np   # Main computations
import numpy.matlib  # np.matlib.repmat - command
import sys           # Error and command line argument handling


#################################################################################
### Parameter definitions
#################################################################################

# "Half-color" in RGB space
h = 128.0/255.0
# Array of available standard colors
colTab = [[1,0,0], # Red
          [1,h,0], # Orange
          [1,1,0], # Yellow
          [0,1,0], # Light Green
          [0,h,0], # Dark Green
          [0,1,1], # Cyan
          [0,0,1], # Blue
          [1,0,1], # Pink
          [h,0,h]] # Purple
colTab = np.array(colTab)
# Dictionary to map the colors to the color indices (internal)
colDict = {'red':1,
           'orange':2,
           'yellow':3,
           'lightgreen':4,
           'darkgreen':5,
           'cyan':6,
           'blue':7,
           'pink':8,
           'purple':9}
# Invert the dictionary to map indices to color-names
colDictInv = {i:c for c,i in colDict.items()}
# rgb_eps (empirical value for python to detect deviations in the colors used and
# still find the "right" color... might need to be adapted)
rgb_eps = 0.5


#################################################################################
### READ IN DATA
#################################################################################

### Get dimensions of vertex and face arrays from .obj file
def getDims(filename):
    vCount=0 # vertex Count = number of vertices
    fCount=0 # face Count = number of faces
    with open(filename,"r") as file:
        # Read v - ...
        line = file.readline(1)
        next(file)
        while(line!='#'):
            vCount = vCount+1
            # Read vn - ...
            line = file.readline(1)
            next(file)
            # Read v - ...
            line = file.readline(1)
            next(file)
        next(file)
        # Read f - ...
        line = file.readline(1)
        while(line):
            next(file)
            fCount = fCount+1
            # Read f - ...
            line = file.readline(1)            
    return vCount, fCount


### Actually read in the vertex and face array
### vertex array contains three information (arrays)
### v = vertex coordinates
### vc = vertex colors
### vn = vertex normals
def getVertAndFaces(filename):
    vCount, fCount = getDims(filename)
    # Vertex array (incl. color info in RGB triple)
    vData = np.zeros((vCount,6),dtype=np.double)
    # Vertex normal array
    vn = np.zeros((vCount,3),dtype=np.double)
    # Connectivity matrix for faces
    t = np.zeros((fCount,3),dtype=np.uint)
    with open(filename,"r") as file:
        for vInd in range(vCount):
            line = file.readline()
            vData[vInd,:] = line.split(' ')[1:]
            line = file.readline()
            vn[vInd,:] = line.split(' ')[1:]
        next(file)
        next(file)
        for fInd in range(fCount):
            line = file.readline()
            line = line.replace('//',' ')
            t[fInd,:] = line.split(' ')[slice(1,6,2)]
        t = t-1
    # Split vertex coordinates from colors    
    v  = np.copy(vData[:,slice(0,3,1)])
    vc = np.copy(vData[:,slice(3,6,1)])
    return v, vc, vn, t


#################################################################################
### MANAGE COLORINGS
#################################################################################

def recolor(vc, colTab, rgb_eps):
    vc_ind = np.zeros((vc.shape[0],1),dtype=np.uint)
    for vInd in range(vc.shape[0]):
        col = vc[vInd,:]
        if np.sum(col)<2.9999999:
            ref_col, ref_col_ind = findColor(col, colTab, rgb_eps)
            vc[vInd,:] = ref_col
            vc_ind[vInd] = ref_col_ind + 1
        else:
            vc[vInd,:] = np.array([1,1,1])
    colDetected = np.setdiff1d(np.unique(vc_ind),np.array([0]))
    print('Colors detected:\n')
    for i in range(colDetected.shape[0]):
        print(str(colDictInv[colDetected[i]]))
    return vc, vc_ind


### Searches the array of available colors <colTab>, which color is closest to <col>
### and gives back that value. If no color is "close enough" it throws an error.
### Close enough is defined by rgb_eps
def findColor(col, colTab, rgb_eps):
    colRep = np.matlib.repmat(col, colTab.shape[0], 1)
    diff = colRep-colTab
    diff_norm = np.linalg.norm(diff, axis=1)
    min_value = np.amin(diff_norm)
    if min_value<rgb_eps:
        ref_col_ind = np.where(diff_norm==min_value)[0]
        ref_col = colTab[ref_col_ind,:]        
    else:
        sys.exit('ERROR: The color '+str(col)+' is ambiguous and can not be matched any default color.')
    return ref_col, ref_col_ind


def colorFaces(t, vc_ind):
    tc_ind = np.zeros((t.shape[0],1),dtype=np.uint)
    for fInd in range(t.shape[0]):
        ind_p = t[fInd,0]
        ind_q = t[fInd,1]
        ind_r = t[fInd,2]
        # Transitive == relation to check if all three vertices have same color
        if (vc_ind[ind_p]==vc_ind[ind_q]) and (vc_ind[ind_q]==vc_ind[ind_r]):
            tc_ind[fInd] = vc_ind[ind_p]
    return tc_ind

#################################################################################
### COMPUTING NEW DATA
#################################################################################

### Computes the face normals via cross product
def getFaceNormals(v, t):
    tn = np.zeros((t.shape[0],3),dtype=np.double)
    for fInd in range(t.shape[0]):
        ind_p = t[fInd,0]
        ind_q = t[fInd,1]
        ind_r = t[fInd,2]
        u = v[ind_q,:]-v[ind_p,:]
        w = v[ind_r,:]-v[ind_p,:]
        n = np.cross(u,w)
        tn[fInd,:] = n/np.linalg.norm(n)
    return tn


#################################################################################
### IDENTIFYING INFLOW (CIRCULAR) BOUNDARY
#################################################################################

### Find (arbitrary) face that is contained in given circle by specifying circle color 
### and searching for a face in that color (which means that all its vertices are of that color)
def findCircFace(tc_ind, colInd):
    for fInd in range(tc_ind.shape[0]):
        if tc_ind[fInd] == colInd:
            return fInd
              

### Create a list of circle boundary points by specifying the circle color and searching
### for all vertices with this color, but within a triangle that also has white vertices
def findCircBoundaryVertList(vc_ind, t, colInd):
    # Initialize with empty length as we do not yet know how long it will become
    cbvl = np.zeros(0, dtype=np.uint)
    for fInd in range(t.shape[0]):
        col_p = vc_ind[t[fInd,0]]
        col_q = vc_ind[t[fInd,1]]
        col_r = vc_ind[t[fInd,2]]
        check = np.array([col_p, col_q, col_r]) == colInd
        if sum(check)==1 or sum(check) == 2:
            if col_p==colInd:
                cbvl=np.append(cbvl, t[fInd,0])
            if col_q==colInd:
                cbvl=np.append(cbvl, t[fInd,1])
            if col_r==colInd:
                cbvl=np.append(cbvl, t[fInd,2])
    # Remove duplicates            
    cbvl = np.unique(cbvl)
    return cbvl


### Select three points from the circle boundary vertex list that are then
### used to compute the circle midpoint and radius
### For now, we just select the first three points in the array. Some more
### elaborate selection strategy might be appropriate
def selectCircleBoundaryPoints(cbvl, v):
    p = v[cbvl[0],:]
    q = v[cbvl[1],:]
    r = v[cbvl[2],:]
    return p, q, r


### Compute the circle midpoint and radius starting from three points on the circle p,q,r and
### the circle plane normal n
def getCircleParams(p,q,r,n):
    a = q-p
    b = r-q
    #c = p-r # Actually not needed
    ma = 1./2*(p+q)
    mb = 1./2*(q+r)
    #mc = 1./2*(r+p) # Actually not needed
    na = np.cross(n,a)
    nb = np.cross(n,b)
    AA = np.column_stack((na,-nb))
    bb = -(ma-mb)
    ATA = np.matmul(np.transpose(AA),AA)
    ATb = np.matmul(np.transpose(AA),bb)
    x = np.linalg.solve(ATA,ATb)
    m = ma+x[0]*na
    rad = np.linalg.norm(m-q)
    # al = np.linalg.norm(a)
    # bl = np.linalg.norm(b)
    # cl = np.linalg.norm(c)
    # rad = (al*bl*cl)/np.sqrt(2.*al**2*bl**2+2.*bl**2*cl**2+2.*cl**2*al**2-al**4-bl**4-cl**4)
    return m, rad


#################################################################################
### MODIFY THE GEOMETRY
#################################################################################

### Translate the object by vector vec by moving all nodes according to it
def translate(v, vec):
    v+=vec
    return v


### Rotation Matrix that maps vector a onto b
### https://math.stackexchange.com/questions/180418/calculate-rotation-matrix-to-align-vector-a-to-vector-b-in-3d
def rotationMatrix(a, b):
    a = a/np.linalg.norm(a)
    b = b/np.linalg.norm(b)
    v = np.cross(a,b)
    #s = np.linalg.norm(v) # Actually not needed
    c = np.dot(a,b)
    if (abs(c+1) < 1e-14):
        sys.exit("Vectors are already aligned in opposite direction!")
    elif (abs(c-1) < 1e-14):
        R = np.eye(3)
    else:
        fct = 1./(1.+c)
        V = np.array([[0,-v[2],v[1]],[v[2],0,-v[0]],[-v[1],v[0],0]])
        R = np.eye(3) + V + fct*np.matmul(V,V)
    return R


### Align geometry by specifying an inflow circle (via its color) and to-be-normal vector
### One then rotates the geometry such that the circle (inflow) normal points in the
### direction we want. Before, the geometry is moved such that the circle mid point gets to
### be at (0,0,0)^T.
def alignGeometry(colInd, nVec, tc_ind, tn, vc_ind, t, v, vn):
    # get current normal of the color-specified circle
    ind = findCircFace(tc_ind, colInd)
    n = tn[ind,:]
    # get circle mid point and radius
    cbvl = findCircBoundaryVertList(vc_ind, t, colInd)
    p,q,r = selectCircleBoundaryPoints(cbvl,v)
    m, rad = getCircleParams(p,q,r,n)
    # Translate geometry
    v = translate(v, -m)
    # Get rotation matrix for orientation rotation
    R = rotationMatrix(n, nVec)
    # Apply rotation
    for vInd in range(v.shape[0]):
        v[vInd,:] = np.matmul(R,v[vInd,:])
        vn[vInd,:] = np.matmul(R,vn[vInd,:])
    for fInd in range(tn.shape[0]):
        tn[fInd,:] = np.matmul(R,tn[fInd,:])
    return v, vn, tn


#################################################################################
### DATA OUTPUT
#################################################################################

def writeOBJ(filename, v, vc, vc_ind, vn, t, colTab):
    with open(filename,"w") as file:
        for vInd in range(v.shape[0]):
            col = vc[vInd,:]
            file.write('v '+str(v[vInd,0])+' '+str(v[vInd,1])+' '+str(v[vInd,2])+' '
                           +str(col[0])+' '+str(col[1])+' '+str(col[2])+'\n')
            file.write('vn '+str(vn[vInd,0])+' '+str(vn[vInd,1])+' '+str(vn[vInd,2])+'\n')
        file.write('# mm_gid 0\n')
        file.write('g mmGroup0\n')
        t = t+1
        for fInd in range(t.shape[0]):
            file.write('f '+str(t[fInd,0])+'//'+str(t[fInd,0])+' '
                           +str(t[fInd,1])+'//'+str(t[fInd,1])+' '
                           +str(t[fInd,2])+'//'+str(t[fInd,2])+'\n')
        

#################################################################################
### MAIN
#################################################################################

# Main 
def Main():
    # Input help
    if len(sys.argv)==1:
        print('How to use: python3 preprocessGeometry "filenameIn.obj" "alignation color" "alignVec"')
        print('For example: python3 preprocessGeometry "AN171.obj", "yellow", "[1,0,0]"')
        sys.exit()
        
    # Argument handling
    filenameIn = sys.argv[1]
    alignCol = sys.argv[2]
    alignColInd = colDict[alignCol]
    alignVecIn = sys.argv[3].split(',')
    alignVec = np.zeros(3)
    alignVec[0] = np.double(alignVecIn[0][1:])
    alignVec[1] = np.double(alignVecIn[1])
    alignVec[2] = np.double(alignVecIn[2][:-1])

    # Obtain data of geometry
    v, vc, vn, t = getVertAndFaces(filenameIn)
    vc, vc_ind = recolor(vc, colTab, rgb_eps)
    tc_ind = colorFaces(t, vc_ind)
    tn = getFaceNormals(v, t)

    # Align geometry
    v, vn, tn = alignGeometry(alignColInd, alignVec, tc_ind, tn, vc_ind, t, v, vn)

    # Write output file
    filenameOut = filenameIn[:-4]+'_aligned.obj'
    writeOBJ(filenameOut, v, vc, vc_ind, vn, t, colTab)


################################################################################
### Entry
################################################################################
if __name__ == "__main__":
    Main()        
