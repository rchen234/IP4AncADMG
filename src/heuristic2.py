import array as arr
import re
import sys
import random
#import pdb; pdb.set_trace()

# TODO: deal with cases 1) when district size > 2 2) nodes in the district being parents of other nodes

def dircyc (n , ne, elist):
    ideg = [0 for i in range(n)]
    odeg = [0 for i in range(n)] #outdegree
    oadjn = [0 for i in range(n*n)] #out adjacency list
    mark = [0 for i in range(n)] #processed nodes
    active = arr.array('i', []) #same as stack
    nchild = arr.array('i', []) #number of processed children
    isactive = [0 for i in range(n)] #idicates if node is on the stack

    nproc = 0
    allcyc = []

    for i in range(ne):
        n1 = elist[2*i]
        n2 = elist[2*i+1]
        oadjn[n1*n + odeg[n1]] = n2
        odeg[n1] = odeg[n1] + 1
        ideg[n2] = ideg[n2] + 1

    while nproc < n: 
        # Find an unmarked node to start from, preferably zero indeg node*/
        start = -1
        pstart = -1
        for i in range(n):
            if mark[i] == 0 and pstart == -1:
                pstart = i
                if ideg[i] == 0:
                    start = i
                    break
        if start == -1 and pstart >= 0:
            start = pstart
        
        #print 'start =', start, ', pstart=',  pstart
        if start == -1:
            break

        active.append(start)
        nchild.append(0)
        mark[start] = 1
        isactive[start] = 1
        
        while len(active) > 0:
            cnode = active[len(active)-1]
            nchildproc = nchild[len(nchild)-1]
            #print 'cnode=',cnode,',nchildproc=',nchildproc, 'act=',len(active), len(nchild)

            if nchildproc == odeg[cnode]:
                nproc = nproc + 1
                isactive[cnode] = 0
                active.pop()
                nchild.pop()
                if len(nchild) > 0:
                    nchild[len(nchild)-1] = nchild[len(nchild)-1] + 1

            else:
                n1 = oadjn[cnode*n + nchildproc]
                #print n1,'mark=',mark[n1]
                if mark[n1] == 0:
                    active.append(n1)
                    nchild.append(0)
                    mark[n1] = 1
                    isactive[n1] = 1
                else:
                    #if n1 in active:
                    if isactive[n1] != 0:
                        # Found previously marked node trace back to cycle
                        pos = len(active)-1
                        cur = active[pos]
                        cyc = [cnode]

                        while cur != n1 and pos > 0:
                            pos = pos -1
                            cur = active[pos]
                            cyc.append(cur)
                        if (cur == n1):
                            # print 'cycle=',cyc
                            allcyc.append(cyc[::-1])
                        else:
                            print('error in computing cycle')
                            exit(-1)

                    # at this point we have a marked child
                    nchild[len(nchild)-1] = nchild[len(nchild)-1] + 1
                    
    return allcyc

def almostdircyc (n , ne, elist, bnode, enode):
    ideg = [0 for i in range(n)]
    odeg = [0 for i in range(n)] #outdegree
    oadjn = [0 for i in range(n*n)] #out adjacency list
    mark = [0 for i in range(n)] #processed nodes
    active = arr.array('i', []) #same as stack
    nchild = arr.array('i', []) #number of processed children
    isactive = [0 for i in range(n)] #idicates if node is on the stack

    nproc = 0
    allcyc = []

    for i in range(ne):
        n1 = elist[2*i]
        n2 = elist[2*i+1]
        oadjn[n1*n + odeg[n1]] = n2
        odeg[n1] = odeg[n1] + 1
        ideg[n2] = ideg[n2] + 1


    mark[bnode] = 1
    active.append(enode)
    nchild.append(0)
    mark[enode] = 1
    isactive[enode] = 1
                
    while len(active) > 0:
        cnode = active[len(active)-1]
        nchildproc = nchild[len(nchild)-1]
        #print 'cnode=',cnode,',nchildproc=',nchildproc, 'act=',len(active), len(nchild)

        if nchildproc == odeg[cnode]:
            nproc = nproc + 1
            isactive[cnode] = 0
            active.pop()
            nchild.pop()
            if len(nchild) > 0:
                nchild[len(nchild)-1] = nchild[len(nchild)-1] + 1

        else:
            n1 = oadjn[cnode*n + nchildproc]
            #print n1,'mark=',mark[n1]
            if mark[n1] == 0:
                active.append(n1)
                nchild.append(0)
                mark[n1] = 1
                isactive[n1] = 1
            else:
                if  n1 == bnode:
                    # Found previously marked node trace back to cycle
                    pos = len(active)-1
                    cur = active[pos]
                    cyc = [cnode]

                    while cur != enode and pos > 0:
                        pos = pos -1
                        cur = active[pos]
                        cyc.append(cur)
                    if (cur == enode):
                        #print 'cycle=',cyc
                        cyc.reverse()
                        cyc.append(bnode)
                        allcyc.append(cyc)
                    else:
                        print('error in computing cycle')
                        exit(-1)
                        
                # at this point we have a marked child
                nchild[len(nchild)-1] = nchild[len(nchild)-1] + 1
                    
    return allcyc

def parse(filename, n, cnodes, cparents, cweight):
    valline = re.compile(':')
    val = re.compile(': ([0-9.]+)')
    parent = re.compile('\(\(([\d+, ]*)\)')
    #parent = re.compile('\(\((\d+[, \d+]*)\)')
    parent2 = re.compile('\(\(([\d+, ]*)\), \(([\d+, ]*)\)')
    #parent2 = re.compile('\(\((\d+[, \d+]*)\), \(([\d+, ]*)\)')
    file = open(filename, "r")
    a = file.readlines()
    
    for l in a:
        if valline.search(l) == None:
            m = parent.match(l).group(1)
            ccomp = [int(num) for num in m.strip(',').split(',')]
        if valline.search(l) != None:
            #print curpar
            ll = l.split(': ')
            pset = [n]
            pset2= [n]
            if parent.match(ll[0]) != None:
                mm = parent.match(ll[0]).group(1)
                if len(mm) != 0:
                    pset = [int(num) for num in mm.strip(',').split(',')]
                if len(ccomp) > 1:
                    mm2 = parent2.match(ll[0]).group(2)
                    #print mm2, len(mm2)
                    if len(mm2) != 0:
                        pset2 = [int(num) for num in mm2.strip(',').split(',')]
                        if len(pset2) == 0:
                            pset2 = [n]
            cweight.append(float(ll[1].strip()))
            cnodes.append(ccomp[:])
            if len(ccomp) > 1:
                cparents.append([pset, pset2])
            else:
                cparents.append([pset])

    #print cnodes
    #print cparents
    #print cweight
    #return cnodes, cparents, cweight

def find_edge(n, wt, elist):
    pos = 0
    e = len(wt)-1
    if e < 0:
        return e
    
    sum = 0.0
    for i in range(len(wt)):
        sum = sum + wt[i]
        
    randval = (sum-1e-6)*random.random()
    #print "r =", randval, sum, e
    
    sum = 0.0
    for i in range(len(wt)):
        sum = sum + wt[i]
        if sum >= randval:
            e = i
            break

    te = e
    while elist[2*te] == n-1 or elist[2*te+1] == n-1:
        te = te - 1
        if te < 0:
            break
    if te < 0:
        te = e
        while elist[2*te] == n-1 or elist[2*te+1] == n-1:
            te = te + 1
            if te >= len(wt):
                break
    return te

def create_undir_edges_from_par(cnodes, cparents, cweight, elist, bdir, wt, inwt):
    edgemap = {}
    for i in range(len(inwt)):
        inwt[i] = 0.0

    for i in range(len(cnodes)):
        if cweight[i] <= 0.001:
            continue
        
        cset = cnodes[i]
        parent = cparents[i]

        if len(cset) > 1:
            for j in cset:
                bdir.append(j)
                
        for j in range(len(cset)):
            node = cset[j]
            par = parent[j]
            if node in par:
                print('mayday', cweight[i], i, cset, par)
                return []
                exit(-1)
            inwt[node] = inwt[node] + cweight[i]
            
            for k in par:
                n1 = k
                n2 = node
                
                if k > node:
                    n1 = node
                    n2 = k
                    
                nodepair = n1*1000 + n2
                if len(edgemap) == 0 or nodepair not in edgemap:
                    edgemap[nodepair] = cweight[i]
                else:
                    edgemap[nodepair] = edgemap[nodepair] + cweight[i]

    #print len(edgemap), edgemap
    for i in edgemap:
        wt.append(edgemap[i])
        n1 = int(i / 1000)
        n2 = int(i % 1000)
        #print n1, n2
        elist.append(n1)
        elist.append(n2)

    #exit(-1)
    #for i in range(len(wt)):
    #    n1 = elist[2*i]
    #    n2 = elist[2*i+1]
    #    print n1, n2, wt[i]

def create_edges_from_par(cnodes, cparents, cweight, elist, bdir, bdirwt, wt, inwt):
    edgemap = {}
    bdirmap = {}
    
    for i in range(len(inwt)):
        inwt[i] = 0.0

    for i in range(len(cnodes)):
        if cweight[i] <= 0.001:
            continue
        
        cset = cnodes[i]
        parent = cparents[i]

        if len(cset) > 1:
            n1 = cset[0]
            n2 = cset[1]
            if n1 > n2:
                temp = n2
                n2 = n1
                n1 = temp
                    
            nodepair = n1*1000 + n2
            if nodepair not in bdirmap:
                bdirmap[nodepair] = cweight[i]
            else:
                bdirmap[nodepair] = bdirmap[nodepair] + cweight[i]
                
        for j in range(len(cset)):
            node = cset[j]
            par = parent[j]
            if node in par:
                print('mayday', cweight[i], i, cset, par)
                return []
                exit(-1)
            inwt[node] = inwt[node] + cweight[i]
            
            for k in par:
                #elist.append(k)
                #elist.append(node)
                #wt.append(cweight[i])

                nodepair = k*1000 + node 
                if nodepair not in edgemap:
                    edgemap[nodepair] = cweight[i]
                else:
                    edgemap[nodepair] = edgemap[nodepair] + cweight[i]

    for i in bdirmap:
        bdirwt.append(bdirmap[i])
        n1 = int(i / 1000)
        n2 = int(i % 1000)
        #print n1, n2
        bdir.append(n1)
        bdir.append(n2)
        
    for i in edgemap:
        wt.append(edgemap[i])
        n1 = int(i / 1000)
        n2 = int(i % 1000)
        #print n1, n2
        elist.append(n1)
        elist.append(n2)

    #print len(elist), elist
    #exit(-1)
    #for i in range(len(wt)):
    #    n1 = elist[2*i]
    #    n2 = elist[2*i+1]
    #    print n1, n2, wt[i]
    
def contract_edge (e, wt, elist, nodelist, cnodes, cparents, cweight, containsbdir):
    node1 = elist[2*e]
    node2 = elist[2*e+1]
    #print node1, node2

    if node1 == node2:
        print('mayday2')
        exit(-1)
    #print node1, node2, len(cnodes), len(cparents), len(cweight)
    #print ' '
    #print 'cnodes = ', cnodes
    #print 'cparents = ', cparents
    #print cweight
    #print 'elist =', elist
    
    for i in nodelist[node1]:
        nodelist[node2].append(i)
    nodelist[node1] = []
    if containsbdir[node1] == 1:
        containsbdir[node2] = 1

    for i in range(len(cnodes)):
        cset = cnodes[i]
        parent = cparents[i]

#        if len(cset) > 1:
#            #assume for now that cset can have only two nodes
#            if node1 in cset and node2 in cset:
#                cnodes[i] = [node2]
#                cparents[i] = list(set(parent[0]).union(parent[1]))
#                containsbdir[node2] = 1
#                if node1 in cparents[i] or node2 in cparents[i]:
#                    cweight[i] = 0.0
#                continue
#            else:
#                if node2 == cset[0] and node1 in parent[0]:
#                    cnodes[i] = [cset[1]]
#                    cparents[i] = [parent[1]]
#                    if node1 in cparents[i]:
#                        ix = cparents[i].index(node1)
#                        cparents[i][ix] = node2
#                    continue
#                else:
#                    if node2 == cset[1] and node1 in parent[1]:
#                        cnodes[i] = [cset[0]]
#                        cparents[i] = [parent[0]]
#                        if node1 in cparents[i]:
#                            ix = cparents[i].index(node1)
#                            cparents[i][ix] = node2
#                        continue
                        
        deln = [0 for j in range(len(cset))]
        ndel = 0
        for j in range(len(cset)):
            
            node = cset[j]
            par = parent[j]
            if node == node2 and node1 in par:
                deln[j] = 1
                ndel = ndel + 1
                #cweight[i] = 0.0                 
                #print('*',cset, par, cweight[i], i)
            else:
               if node == node1:
                   if node2 in par:
                       deln[j] = 1
                       #cweight[i] = 0.0                
                       #print('+',cset, par)
                   else:
                       #print('bb',cset, par)
                       cset[j] = node2
                       #print('aa',cset, par)
               else:       
                   if node1 in par:
                       ix = par.index(node1)
                       if node2 not in par:
                           #print('b',cset, par, cweight[i], i)
                           par[ix] = node2
                           #print('a',cset, par, cweight[i], i)
                       else:
                           #print('br',cset, par)
                           par.remove(node1)
                           #print('br',cset, par)

        
        if ndel > 0:
            if ndel == len(cset):
                cweight[i] = 0.0
            else:
                newcset = []
                newpar = []
                for j in range(len(cset)):
                    if deln[j] == 0:
                        newcset.append(cset[j])
                        newpar.append(parent[j])
                cnodes[i] = newcset
                cparents[i] = newpar
                del deln
        #end for j
    #end for i
    
def contract_bdir_edge (node1, node2, wt, elist, nodelist, cnodes, cparents, cweight, containsbdir):
    if node1 == node2:
        print('mayday2')
        exit(-1)

    #print node1, node2, len(cnodes), len(cparents), len(cweight)
    #print ' '
    #print 'cnodes = ', cnodes
    #print 'cparents = ', cparents
    #print cweight
    #print 'elist =', elist
    
    for i in nodelist[node1]:
        nodelist[node2].append(i)
    nodelist[node1] = []
    containsbdir[node1] = 0
    containsbdir[node2] = 1

    for i in range(len(cnodes)):
        cset = cnodes[i]
        parent = cparents[i]

        if len(cset) == 1:
            if node1 in cset or node2 in cset:
                cweight[i] = 0.0
                continue
        if len(cset) > 1:
            if node1 in cset and node2 not in cset:
                cweight[i] = 0.0
                continue
            if node2 in cset and node1 not in cset:
                cweight[i] = 0.0
                continue
            if node1 in cset and node2 in cset:
                cnodes[i] = [node2]
                cparents[i] = [list(set(parent[0]).union(parent[1]))]
                continue
                                    
    for i in range(len(cnodes)):
        cset = cnodes[i]
        parent = cparents[i]

        if node1 in cset or node2 in cset:
            continue
        for j in range(len(cset)):
            node = cset[j]
            par = parent[j]
            if node1 in par:
                ix = par.index(node1)
                if node2 not in par:
                    #print('b',cset, par, cweight[i], i)
                    par[ix] = node2
                    #print('a',cset, par, cweight[i], i)
                else:
                    #print('br',cset, par)
                    par.remove(node1)
                    #print('br',cset, par)

def single_contraction_round(n, cnodes, cparents, cweight, allcyc, cutwt, dbg):
    elist = arr.array('i',[]) #array to store edges
    bdir = arr.array('i',[]) #array to store bidirected edges
    wt = arr.array('d',[]) #array to store edges
    inwt = arr.array('d', [0.0 for i in range(n)])
    nodelist = [[i] for i in range(n)]
    containsbdir = [0 for i in range(n)]

    create_undir_edges_from_par(cnodes, cparents, cweight, elist, bdir, wt, inwt)
    #print 'bef', len(cnodes), len(cweight), inwt
    over = n-5
    while over >= 0:
        e = find_edge(n, wt, elist)
        
        if e < 0 or e >= len(wt):
            #print e
            #print ('no edge to contract')
            over = -1
            break
        if dbg > 0:
            print('e =', e, 'len =', len(wt), 'n1,n2=', elist[2*e], elist[2*e+1])
            for i in range(len(cnodes)):
                print(cnodes[i], cparents[i], cweight[i])
        contract_edge (e, wt, elist, nodelist, cnodes, cparents, cweight, containsbdir)
        
        del wt
        del elist
        del bdir
        wt = arr.array('d',[]) #array to store edges
        elist = arr.array('i',[]) #array to store edges
        bdir = arr.array('i',[]) #array to store bidirected edges

        create_undir_edges_from_par(cnodes, cparents, cweight, elist, bdir, wt, inwt)
        
        for i in range(len(inwt)):
            if inwt[i] < 0.98*cutwt and len(nodelist[i]) > 1 and nodelist[i] not in allcyc:
                allcyc.append(nodelist[i][:])
                #print('hello', nodelist[i],  inwt[i])

        if dbg > 0:        
            print (nodelist)
            print (inwt)
        over = over - 1
    #end while
#end function

def contract_heur(n, icnodes, icparents, icweight, cutwt):
    #first locate edges with weight>= 0.99 and check if almost directed cycle exists
    wt = arr.array('d',[]) #array to store edges
    elist = arr.array('i',[]) #array to store edges
    bdir = arr.array('i',[]) #array to store bidirected edges
    bdirwt = arr.array('d',[]) #array to store bidirected edges
    #telist = arr.array('i',[]) #array to store edges
    allcyc = []
    inwt = arr.array('d', [0 for i in range(n)])
    cnodes = []
    cparents = []
    cweight = []
    
    # copy all input arrays as we will modify them
    for i in range(len(icnodes)):
        if icweight[i] > 0.001:
            cnodes.append(icnodes[i][:])
            cparents.append(icparents[i][:])
            cweight.append(icweight[i])
    
    nccomp = len(cnodes)
    #print('n = ', n, 'nccomp = ', nccomp)

    #for k in range(len(cnodes)):
    #    print 'he', k, cnodes[k], cparents[k], cweight[k]

    create_edges_from_par(cnodes, cparents, cweight, elist, bdir, bdirwt, wt, inwt)
    
    '''
    #now create array of wt 0.95+ edges and test for directed cycles
    for i in range(len(wt)):
        if wt[i] >= 0.95:
            telist.append(elist[2*i])
            telist.append(elist[2*i+1])

    allcyc = dircyc(n, len(telist)/2, telist)
    if len(allcyc) > 0:
        return allcyc
    # at this point we do not have any directed cycles
    '''

    '''
    # now look for cluster ineqs. First modify bidirected sets, here we assume we are splitting it up
    if len(bdir) > 0:
        cl = len(cnodes)
        for i in range(cl):
            cset = cnodes[i]
            parent = cparents[i]
            twt = cweight[i]

            if len(cset) > 1:
                #print (cset, len(cset))
                cnodes[i] = [cset[0]]
                cparents[i] = [parent[0][:]]
                cnodes.append([cset[1]])
                cweight.append(twt)
                cparents.append([parent[1][:]])
        #print ('old len = ', cl, ' new len = ', len(cnodes))

        i = len(cnodes)-1
        while i >=0:
            if cweight[i] <= 0.001:
                del cnodes[i]
                del cparents[i]
                del cweight[i]
            i = i-1
    '''
    
    for numrep in range(20):
        dbg = 0
        # now modify cnodes, cparents, cweight to find regular cluster inequalities
        tnodes = []
        tparents = []
        tweight = arr.array('d', [])

        # now copy all parent set info as it is used destructively inside single_contraction round
        for i in range(len(cnodes)):
            cset = []
            par = []
            for j in range(len(cnodes[i])):
                cset.append(cnodes[i][j])

            for j in range(len(cparents[i])):
                par.append(cparents[i][j][:])

            tnodes.append(cset)
            tparents.append(par)
            tweight.append(cweight[i])

        if dbg > 0:
            for i in range(len(tnodes)):
                print(cnodes[i], cparents[i], cweight[i])

        single_contraction_round (n, tnodes, tparents, tweight, allcyc, cutwt, dbg)
        
        if len(allcyc) > 10:
            break
    #end repeat
    
    return allcyc
#end of function

def contract_heur_bdir(n, icnodes, icparents, icweight):
    #first locate edges with weight>= 0.99 and check if almost directed cycle exists
    wt = arr.array('d',[]) #array to store edges
    elist = arr.array('i',[]) #array to store edges
    bdir = arr.array('i',[]) #array to store bidirected edges
    bdirwt = arr.array('d',[]) #array to store bidirected edges
    inwt = arr.array('d', [0 for i in range(n)])
    cnodes = []
    cparents = []
    cweight = []
    allcyc = []
    
    # copy all input arrays as we will modify them
    for i in range(len(icnodes)):
        if icweight[i] > 0.001:
            cnodes.append(icnodes[i][:])
            cparents.append(icparents[i][:])
            cweight.append(icweight[i])
    
    nccomp = len(cnodes)
    #print('n = ', n, 'nccomp = ', nccomp)

    
    create_edges_from_par(cnodes, cparents, cweight, elist, bdir, bdirwt, wt, inwt)
    
    for b in range(int(len(bdir)/2)):

        # now modify cnodes, cparents, cweight to find regular cluster inequalities
        tnodes = []
        tparents = []
        tweight = arr.array('d', [])
        nodelist = [[i] for i in range(n)]
        containsbdir = [0 for i in range(n)]

        #print('bdir = ', bdir[2*b], bdir[2*b+1], bdirwt[b])

        for i in range(len(cnodes)):
            cset = []
            par = []
            for j in range(len(cnodes[i])):
                cset.append(cnodes[i][j])

            for j in range(len(cparents[i])):
                par.append(cparents[i][j][:])

            tnodes.append(cset)
            tparents.append(par)
            tweight.append(cweight[i])


        contract_bdir_edge (bdir[2*b], bdir[2*b+1], wt, elist, nodelist, tnodes, tparents, tweight, containsbdir)

        #for k in range(len(tnodes)):
        #    print 'll', k, tnodes[k], tparents[k], tweight[k]

        newcyc = contract_heur (n, tnodes, tparents, tweight, bdirwt[b])
        if len(newcyc) > 0:
            for cyc in newcyc:
                if bdir[2*b+1] in cyc:
                    cyc1 = cyc
                    cyc1.remove(bdir[2*b+1])
                    cyc1.sort()
                    newc = [bdir[2*b], bdir[2*b+1]] + cyc

                    if not newc in allcyc:
                        allcyc.append(newc)
                        #print('adding ', newc)
        
    return allcyc
#end of function

def add_dummy_node(n, cnodes, cparents, cweight):
    for i in range(len(cnodes)):
        cset = cnodes[i]
        parent = cparents[i]

        for j in range(len(cset)):
            if len(parent[j]) == 0:
                parent[j].append(n)
                
#n = 3
#ne = 3
#elist = arr.array('i', [0,1, 1, 2, 2, 0])
#n = 6
#ne = 7
#elist = arr.array('i', [0,1, 1, 2, 2,0, 3,4,4,5,5,3, 3,1])
#allc = almostdircyc(n, ne, elist, 3, 2)
#print(allc)
'''
#parse("tst1.sol", 20)
nco = [7, 6, 11, 6, 0, 6, 13, 6, 10, 6, 15, 6, 4, 6]
ncount = 0

n = 5
gcnodes = [[0], [1], [2], [3,4]]
gcweight = [1.0, 1.0, 1.0, 1.0]
gcparents = [[[1,2]], [[2, 3]], [[]], [[], [1]]]

#add_dummy_node(n, gcnodes, gcparents, gcweight)
#print gcnodes
#print gcparents
             
gcnodes = []
gcparents = []
gcweight = arr.array('d', [])

bdiredgewt = 0.0

if len(sys.argv) > 1:
    n=int(sys.argv[1])
    filename = sys.argv[2]

    parse(filename, n, gcnodes, gcparents, gcweight)

random.seed(1)
rset = contract_heur_bdir(n+1, gcnodes, gcparents, gcweight )
print(rset)
'''