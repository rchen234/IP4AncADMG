import array as arr

def dircyc (n , ne, elist):
    deg = [0 for i in range(n)]
    ideg = [0 for i in range(n)]
    odeg = [0 for i in range(n)]
    nchild = [0 for i in range(n)]
    bdeg = [0 for i in range(n)]
    adjn  = [0 for i in range(n*n)]
    oadjn = [0 for i in range(n*n)]
    iadjn = [0 for i in range(n*n)]
    mark = [0 for i in range(n)]
    process = [0 for i in range(n)]
    curnode = arr.array('i', [])
    prev = [-1 for i in range(n)]
    
    allcyc = []

    nproc = 0
    for i in range(ne):
        n1 = elist[2*i]
        n2 = elist[2*i+1]
        oadjn[n1*n + odeg[n1]] = n2
        iadjn[n2*n + ideg[n2]] = n1
        odeg[n1] += 1
        ideg[n2] += 1
        adjn[n1*n + deg[n1]] = n2
        adjn[n2*n + deg[n2]] = n1
        deg[n1] += 1
        deg[n2] += 1

    while nproc < n: 
        # Find an unmarked node to start from, preferably zero indeg node*/
        start = -1
        pstart = -1
        for i in range(n):
            if process[i] == 0 and pstart == -1:
                pstart = i
                if ideg[i] == 0:
                    start = i
                    break
        if start == -1 and pstart >= 0:
            start = pstart
        
        #print('start =', start, ', pstart=',  pstart)
        if start == -1:
            break

        curnode.append(start)
        mark[start] = 1
        nproc += 1
        
        while len(curnode) > 0:
            cnode = curnode.pop()
            #curnode[len(curnode)-1]
            #print('cnode=',cnode)
            mark[cnode] = 1
            
            for i in range(odeg[cnode]):
                n1 = oadjn[cnode*n + i]
                #print(n1,'mark=',mark[n1])
                if mark[n1] > 0 and process[n1] == 0:
                    # Found previously marked node trace back to cycle
                    cur = cnode
                    cyc = [cnode]
                    while cur != n1 and cur != -1:
                        cur = prev[cur]
                        cyc.append(cur)
                    if (cur == n1):
                        #print('cycle=',cyc)
                        allcyc.append(cyc)
                if mark[n1] == 0 and process[n1] == 0:        
                    curnode.append(n1)
                    prev[n1] = cnode
                    
        for i in range(n):
            if mark[i] == 1:
                process[i] = 1     

    return allcyc

#n = 3
#ne = 3
#elist = arr.array('i', [0,1, 1, 2, 2, 0])
'''
n = 6
ne = 7
elist = arr.array('i', [0,1, 1, 2, 0,2, 3,4,4,5,5,3, 3,1])
allc = dircyc(n, ne, elist)
print(allc)
'''