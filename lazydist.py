
import numpy as np

def aargmax(x):
    N = np.argmax(x)
    i = N//x.shape[1]
    j = N%x.shape[1]
    return x[i,j],i,j

def dstore(x,i,j,v):
    # Store newly calculated distance v between items i,j in matrix x
    x[i,j] = v
    x[j,i] = v

def propagateLO(prob,i,j):
    # The lower bound for d(x[i],x[j]) has been updated (made higher)
    # i.e. prob.lo[i,j] > prob.lplo[i,j]
    # This has consequences for other bounds that can be propagated
    # (see README)

    dstore(prob.lplo,i,j,prob.lo[i,j])

    #######################################################
    # For every a,b, we want to replace lo[a,b] with      #
    # lo[i,j] - hi[i,a] - hi[j,b]                         #
    # if this is larger than the current value of lo[a,b] #
    #######################################################

    # A is a vector (scalar-vector arithmetic occurs entry-wise)
    A = prob.lo[i,j] - prob.hi[i]
    # B is a matrix with constant columns
    B = np.tile(prob.hi[j],(prob.n,1)).T
    # vector-matrix arithmetic occurs row-wise
    newlo = A - B
    # Third arg of maximum is storage location
    np.maximum(prob.lo,newlo  ,prob.lo)
    np.maximum(prob.lo,newlo.T,prob.lo) # Reverse roles of A and B

def propagateHI(prob,i,j):
    # The upper bound for d(x[i],x[j]) has been updated (made lower)
    # i.e. prob.hi[i,j] < prob.lphi[i,j]
    # This has consequences for other bounds that can be propagated
    # (see README)

    dstore(prob.lphi,i,j,prob.hi[i,j])

    ########################################################
    # For every a,b, we want to replace hi[a,b] with       #
    # hi[a,i] + hi[i,j] + hi[j,b]                          #
    # if this is smaller than the current value of hi[a,b] #
    ########################################################

    # See line-by-line comments in propagateLO
    A = prob.hi[i,j] + prob.hi[i]
    B = np.tile(prob.hi[j],(prob.n,1)).T
    newhi = A + B
    np.minimum(prob.hi,newhi  ,prob.hi)
    np.minimum(prob.hi,newhi.T,prob.hi)

    #######################################################
    # and we want to replace lo[a,b] with                 #
    # lo[i,a] - hi[i,j] - hi[j,b]                         #
    # if this is larger than the current value of lo[a,b] #
    #######################################################

    A = prob.lo[i] - prob.hi[i,j]
    B = np.tile(prob.hi[j],(prob.n,1)).T
    newlo = A - B
    np.maximum(prob.lo,newlo  ,prob.lo)
    np.maximum(prob.lo,newlo.T,prob.lo)

def fillin(prob):
    # Measure some distances to tighten upper and lower bounds
    # until tolerance is met

     ###########################################
     # Measure from everything to a random hub #
     ###########################################
    hub = np.random.randint(prob.n)
    for i in xrange(prob.n):
        if i == hub:
            continue
        dhub = prob.d(prob.x[hub],prob.x[i])
        dstore(prob.hi,hub,i,dhub)
        dstore(prob.lo,hub,i,dhub)

    # Manually set the upper and lower bounds from the hub distances
    dd = np.tile(prob.hi[hub],(prob.n,1))
    prob.lo = abs(dd-dd.T)
    prob.hi = dd+dd.T
    prob.hi[xrange(prob.n),xrange(prob.n)] = 0

    while True:
        ######################
        # Measure a distance #
        ######################
        c = prob.crit(prob.lo,prob.hi)
        if prob.toltype == 'local':
            c[xrange(prob.n),xrange(prob.n)] = 0
        cij,i,j = aargmax(c)
#        print cij,i,j
        # If max criteria is below tol, return
        if cij < prob.tol:
            return
        # Otherwise, measure distance where criteria is greatest
        # (which reduces uncertainty to 0)
        dij = prob.d(prob.x[i],prob.x[j])
        dstore(prob.hi,i,j,dij)
        dstore(prob.lo,i,j,dij)
        while True:
            #########################
            # Propagate information #
            #########################
            # Choose bound that has been updated most since last propagation
            dhi,ihi,jhi = aargmax(prob.lphi-prob.hi  )
            dlo,ilo,jlo = aargmax(prob.lo  -prob.lplo)
            # If everything is up to date, break
            if dhi == 0 and dlo == 0:
                break
            # Otherwise, propagate most out-of-date information
            if dhi > dlo:
                propagateHI(prob,ihi,jhi)
            else:
                propagateLO(prob,ilo,jlo)
class LazyDistProb:
    # Having a thin wrapper class reduces the length of parameter lists
    def __init__(self,x,d,tol,toltype):
        if not toltype in ['global','local']:
            raise ValueError('Unknown tolerance type: {}'.format(toltype))

        # Lower bounds
        lo = np.zeros([n,n],dtype=float)
        # Upper bounds
        hi = np.empty([n,n],dtype=float)
        hi[:,:] = np.inf
        hi[xrange(n),xrange(n)] = 0
        if toltype == 'local':
            def crit(lo,hi):
                # Upper-lower bound gap
                # relative to upper bound
                return (hi-lo)/hi
        else:
            def crit(lo,hi):
                # Upper-lower bound gap
                # relative to lower bound for longest distance
                return (hi-lo)/(lo.max())
        self.crit = crit

        self.n = len(x)
        self.x=x
        self.d=d
        self.lo=lo
        self.hi=hi
        self.tol=tol
        self.toltype=toltype
        # value of upper bounds the last time a change was propagated
        self.lphi=hi.copy()
        # same but for lower bounds
        self.lplo=lo.copy()

def lazyDist(x,d,lo=None,hi=None,tol=1e-3,toltype='global'):
    '''
    x: integer-indexed iterable of objects to compare
    d: a metric function
        d(x[i],x[j]) >= 0 should act like a proper metric
    lo,hi: current upper and lower bounds (0 and inf if no others can be provided)
        if the distance i,j has already been calculated,
        set lo[i,j]=lo[j,i]=hi[i,j]=hi[j,i] = distance i,j
    tol: tolerance
    toltype:
    'global' each distance is calculated within a proportion of longest distance
    'local' each distance is calculated within a proportion of its upper bound
    '''
    if not lo is None or not hi is None:
        raise ValueError('Sorry, a priori bounds not implemented yet.')
    prob = LazyDistProb(x,d,tol,toltype)
    fillin(prob)
    return prob.lo,prob.hi

if __name__=='__main__':
    np.random.seed(16)
    ngroups = 4
    npergroup = 5
    n = ngroups*npergroup
    x = np.random.randn(n) + 1j*np.random.randn(n)
    x[  npergroup:2*npergroup] += 20
    x[2*npergroup:3*npergroup] += 20j
    x[3*npergroup:4*npergroup] += 20+20j

    class DistCounter:
        def __init__(self):
            self.calcs = []
        def __call__(self,a,b):
            self.calcs.append( (a,b) )
            return abs(a-b)

    d = DistCounter()

    tol = 1e-1
    lo,hi = lazyDist(x,d,tol=tol,toltype='local')
    xx = np.tile(x,(n,1))
    dd = abs(xx-xx.T)
    print 'Calculated {} out of {} distances'.format(len(d.calcs),n*(n-1)/2)
    print 'Checking that actual distances are in correct range'
    print (dd>=lo-1e6).all()
    print (dd<=hi+1e6).all()
    cc = (hi-lo)/hi
    cc[xrange(n),xrange(n)] = 0
    print (cc<tol).all()

    import matplotlib.pyplot as plt

    for a,b in d.calcs:
        plt.plot([a.real,b.real],[a.imag,b.imag],'ko-')
    plt.show()

