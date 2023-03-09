import math
import random

the = {
    'bootstrap':512,
    'conf':0.05,
    'cliff':.4,
    'cohen':.35,
    'Fmt':'%6.2f',
    'width':40
}


def erf(x):
    # from Abramowitz and Stegun 7.1.26 
    # https://s3.amazonaws.com/nrbook.com/AandS-a4-v1-2.pdf
    # (easier to read at https://en.wikipedia.org/wiki/Error_function#Approximation_with_elementary_functions)
    a1 =  0.254829592
    a2 = -0.284496736
    a3 =  1.421413741
    a4 = -1.453152027
    a5 =  1.061405429
    p  =  0.3275911
    # Save the sign of x
    sign = 1
    if x < 0:
        sign = -1
    x = abs(x)
    # A&S formula 7.1.26
    t = 1.0/(1.0+p*x)
    y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*math.exp(-x*x)
    return sign*y

def gaussian(mu=0, sd=1):
    # --> n; return a sample from a Gaussian with mean `mu` and sd `sd`
    return mu + sd * math.sqrt(-2*math.log(random.random())) * math.cos(2*math.pi*random.random())

def samples(t, n=None):
    u = {}
    if len(t)!=0:
        for i in range(1, n+1 if n is not None else len(t)+1):
            u[i]=int(len(t)*random.random()) # this is fast but not accurate
            # u[i]=random.choice(list(t.values())) # this is extreme slow but gives correct answer
    return u

def cliffsDelta(ns1, ns2):
    n = 0
    gt = 0
    lt = 0
    if len(ns1) > 128:
        ns1 = samples(ns1, 128)
    if len(ns2) > 128:
        ns2 = samples(ns2, 128)
    for _,x in ns1.items():
        for _,y in ns2.items():
            n += 1
            if x > y: gt += 1
            if x < y: lt += 1
    return abs(lt - gt)/n <= the['cliff']

def add(i, x):
    i['n'] += 1
    d = x-i['mu']
    i['mu'] += d/i['n']
    i['m2'] += d*(x-i['mu'])
    i['sd'] = 0 if i['n']<2 else math.sqrt(i['m2']/(i['n']-1))

def NUM(t=None):
    i = {
        'n':0,
        'mu':0,
        'm2':0,
        'sd':0
    }
    for _,x in t.items() if t is not None else {}.items():
        add(i, x)
    return i

def delta(i, other):
    e = 1E-32
    y = i
    z = other
    return abs(y['mu'] - z['mu']) / math.sqrt(e + y['sd']**2/y['n'] + z['sd']**2/z['n'])

def bootstrap(y0, z0):
    x, y, z, yhat, zhat = NUM(), NUM(), NUM(), {}, {}
    # x will hold all of y0,z0
    # y contains just y0
    # z contains just z0
    for _,y1 in y0.items():
        add(x, y1)
        add(y, y1)
    for _,z1 in z0.items():
        add(x, z1)
        add(z, z1)
    xmu, ymu, zmu = x['mu'], y['mu'], z['mu']
    # yhat and zhat are y,z fiddled to have the same mean
    for _,y1 in y0.items():
        yhat[len(yhat)+1] = y1 - ymu + xmu
    for _,z1 in z0.items():
        zhat[len(zhat)+1] = z1 - zmu + xmu
    # tobs is some difference seen in the whole space
    tobs = delta(y, z)
    n = 0
    for _ in range(1, the['bootstrap']+1):
        # here we look at some delta from just part of the space
        # it the part delta is bigger than the whole, then increment n
        if delta(NUM(samples(yhat)), NUM(samples(zhat))) > tobs:
            n += 1
    # if we have seen enough n, then we are the same
    # On Tuesdays and Thursdays I lie awake at night convinced this should be "<"
    # and the above "> obs" should be "abs(delta - tobs) > someCriticalValue". 
    return n / the['bootstrap'] >= the['conf']

def dictsort(d):
    sorteddict = sorted(d.items(), key=lambda item:item[1])
    u = {}
    n = 0
    for k,v in sorteddict:
        n+=1
        u[n] = v
    return u


def RX(t, s):
    dictsort(t)
    return {
        'name':s if s is not None else "",
        'rank':0,
        'n':len(t),
        'show':'',
        'has':t
    }

def mid(t):
    t = t['has'] if 'has' in t.keys() else t
    n = len(t)//2
    return (t[n] + t[n+1])/2 if len(t)%2==0 else t[n+1]

def div(t):
    t = t['has'] if 'has' in t.keys() else t
    return (t[len(t)*9//10] - t[len(t)*1//10])/2.56

def merge(rx1, rx2):
    rx3 = RX({}, rx1['name'])
    for _,t in rx1['has'].items():
        for _,x in t.items():
            rx3['has'][len(rx3['has'])+1] = x
    for _,t in rx2['has'].items():
        for _,x in t.items():
            rx3['has'][len(rx3['has'])+1] = x
    rx3['has'] = dictsort(rx3['has'])
    rx3['n'] = len(rx3['has'])
    return rx3

def scottKnot(rxs):
    def merges(i, j):
        out = RX({},rxs[i]['name'])
        for k in range(i,j+1):
            out = merge(out, rxs[j]) # shouldn't here be k?
        return out
    def same(lo, cut, hi):
        l = merges(lo, cut)
        r = merges(cut+1, hi)
        return cliffsDelta(l['has'],r['has']) and bootstrap(l['has'], r['has'])
    
    def recurse(lo, hi, rank):
        b4 = merges(lo, hi)
        best = 0
        for j in range(lo,hi+1):
            if j < hi:
                l = merges(lo, j)
                r = (j+1, hi)
                now = (l['n']*(mid(l) - mid(b4))**2 + r['n']*(mid(r) - mid(b4))**2) / (l['n'] + r['n'])
                if now > best:
                    if abs(mid(l) - mid(r)) >= cohen:
                        cut, best = j ,now
        if cut and not same(lo,cut,hi):
            rank = recurse(lo, cut, rank) + 1
            rank = recurse(cut+1, hi ,rank)
        else:
            for i in range(lo, hi+1):
                rxs[i]['rank']=rank
    sortedrxs = sorted(rxs.items(), key = lambda item:mid(item[1]))
    newrxs = {}
    n=0
    for k,v in sortedrxs:
        n+=1
        newrxs[n]=v
    cohen = div(merges(1, len(newrxs))) * the['cohen']
    recurse(1, len(rxs), 1)
    return newrxs

def tiles(rxs):
    lo, hi = -math.inf, math.inf
    for _,rx in rxs.items():
        lo = math.min(lo,rx['has'][1])
        hi = math.hi(hi, rx['has'][len(rx['has'])])
    for _,rx in rxs.items():
        t = rx['has']
        u = {}
        def of(x, most):
            return math.max(1, math.min(most, x))
        def at(x):
            return t[of(int(len(t)*x), len(t))]
        def pos(x):
            return math.floor(of(the['width']*(x-lo)/int(hi-lo+1E-32), the['width']))
        for i in range(1, the['width']+1):
            u[1+len(u)]=" "
        a,b,c,d,e = at(.1), at(.3), at(.5), at(.7), at(.9) 
        A,B,C,D,E = pos(a), pos(b), pos(c), pos(d), pos(e)
        for i in range(A, B+1):
            u[i]="-"
        for i in range(D, E+1):
            u[i]="-"
        u[the['width']//2] = "|"
        u[C] = "*"
        import string
        rx['show'] = str(u)+"{"+string.format(the['fmt'],a)
        for _,x in {b,c,d,e}.items():
            rx['show'] += "," + string.format(the['fmt'],x)
        rx['show'] += "}"
    return rxs




eg = {}


def okfun(n=1):
    random.seed(n if n is not None else 1)
eg['ok'] = okfun

def numfun():
    n = NUM({1:1,2:2,3:3,4:4,5:5,6:6,7:7,8:8,9:9,10:10})
    print("  {}  {}  {}".format(n['n'],n['mu'],n['sd']))
eg['num'] = numfun

def gaussfun():
    t = {}
    for i in range(1,10**4+1):
        t[1+len(t)]=gaussian(10,2)
    n = NUM(t)
    print("  {}  {}  {}".format(n['n'],n['mu'],n['sd']))
eg['gauss'] = gaussfun

def bootmufun():
    a = {}
    for i in range(1,10**4+1):
        a[1+len(a)]=gaussian(10,1)
    print("  "+"mu "+"sd "+"cliffs "+"boot "+"both ")
    print("  "+"-- "+"-- "+"------ "+"---- "+"---- ")
    for mu in range(100,111):
        mu=float(format(mu*0.1, '.1f'))
        b={}
        for i in range(1,101):
            b[1+len(b)]=gaussian(mu,1)
        cl = cliffsDelta(a, b)
        bs = bootstrap(a, b)
        print(" {} 1 {}  {} {}".format(mu,cl,bs,cl and bs))
# eg['bootmu'] = bootmufun

def basicfun():
    print("\t\tTrue \t{}\t{}".format(bootstrap({1:8,2:7,3:6,4:2,5:5,6:8,7:7,8:3},
                                            {1:8,2:7,3:6,4:2,5:5,6:8,7:7,8:3}),
                                cliffsDelta({1:8,2:7,3:6,4:2,5:5,6:8,7:7,8:3},
                                            {1:8,2:7,3:6,4:2,5:5,6:8,7:7,8:3})))
    print("\t\tTrue \t{}\t{}".format(bootstrap({1:8,2:7,3:6,4:2,5:5,6:8,7:7,8:3},
                                            {1:9,2:9,3:7,4:8,5:10,6:9,7:6}),
                                cliffsDelta({1:8,2:7,3:6,4:2,5:5,6:8,7:7,8:3},
                                            {1:9,2:9,3:7,4:8,5:10,6:9,7:6})))
    print("\t\tFalse\t{}\t{}".format(bootstrap({1:0.34, 2:0.49, 3:0.51, 4:0.6,   5:.34,  6:.49,  7:.51, 8:.6}, 
                                             {1:0.6,  2:0.7,  3:0.8,  4:0.9,   5:.6,   6:.7,   7:.8,  8:.9}),
                                 cliffsDelta({1:0.34, 2:0.49, 3:0.51, 4:0.6,   5:.34,  6:.49,  7:.51, 8:.6}, 
                                             {1:0.6,  2:0.7,  3:0.8,  4:0.9,   5:.6,   6:.7,   7:.8,  8:.9})))
eg['basic'] = basicfun

def prefun():
    print('\neg3')
    d = 1
    for i in range(1,11):
        t1,t2={},{}
        for j in range(1,33):
            t1[1+len(t1)]=gaussian(10,1)
            t2[1+len(t2)]=gaussian(d*10,1)
        d=float(format(d,'.2f'))
        print("\t{}\t{} {} {}".format(d,d<1.1,bootstrap(t1,t2),bootstrap(t1,t1)))
        d+=0.05
        
eg['pre']=prefun

for k, fun in eg.items():
    eg['ok']()
    print("\n"+str(k))
    fun()
