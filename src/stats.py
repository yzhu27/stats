import math


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
    return mu + sd * math.sqrt(-2*math.log(math.random())) * math.cos(2*math.pi*math.random())

def samples(t, n):
    u = {}
    for i in range(1, n+1 if n is not None else len(t+1)):
        u[i]=t[math.random(len(t))]
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

def NUM(t):
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
    return math.abs(y['mu'] - z['mu']) / math.sqrt(e + y['sd']**2/y['n'] + z['sd']**2/z['n'])

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
    sorteddict = sorted(d.items(), key=lambda item:item[0])
    u = {}
    for k,v in sorteddict:
        u[k] = v
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

    