#this function return the value of the input formula("a") depending on our order("order") as boolean
def answer(a,binarylist,ordered):
    if(len(a)==1):
        index1=ordered.index(a[0])
        if(binarylist[index1]==0):
            return True
        else:
            return False
    else:
        if(a[1]=='~'):
            par=0
            for i in range(1,len(a)-1):
                if(a[i]=='('):
                    par+=1
                elif(a[i]==')'):
                    par-=1
            if(par==0):
                return (not answer(a[2:len(a)-1],binarylist,ordered))
        else:
            par=0
            index=0
            for i in range(1,len(a)-1):
                if(a[i]=='('):
                    par+=1
                elif(a[i]==')'):
                    par-=1
                if (par==0 and (a[i]=='&' or a[i]=='|' or a[i]=='>')):
                    index = i
                    break
            if(a[index]=='&'):
                return (answer(a[1:index],binarylist,ordered) and answer(a[index+1:len(a)-1],binarylist,ordered))
            elif(a[index]=='|'):
                return (answer(a[1:index], binarylist,ordered) or answer(a[index + 1:len(a) - 1], binarylist,ordered))
            elif(a[index]=='>'):
                return (not answer(a[1:index], binarylist,ordered) or answer(a[index + 1:len(a) - 1], binarylist,ordered))
        return False
#end

#this function increase a list like binary number (example before:"[0,1,0,1]" after:"[0,1,1,0])
def generate(binarylist):
    binarylist[-1]+=1
    for i in range((len(binarylist)-1),0,-1):
        if(binarylist[i]>1):
            binarylist[i]=0
            binarylist[i-1]+=1
#end

#this function generate the truth table of all of the possible interpolators.
#return a list of truth table of all of the possible interpolators.
#we have three case for each row
#case1: if A=1 and B=1 so C=1
#case2: if A=0 and B=0 so C=0
#case3: if A=0 and B=1 so C=0 or C=1 (C can get the both value)
def generate_inter(truthtable_a,truthtable_b):
    interpolator=[""]
    for i in range(len(truthtable_a)):
        if(truthtable_a[i]==truthtable_b[i]):
            for j in range(len(interpolator)):
                interpolator[j]+=truthtable_a[i]
        else:
            tmp=interpolator.copy()
            for j in range(len(interpolator)):
                interpolator[j]+="F"
            for j in range(len(interpolator)):
                tmp[j]+="T"
            interpolator=interpolator+tmp
    return interpolator
#end

#this function get the one row of truth table that is 1 and make the minterm of that row
def binary_to_minterm(binaryvalue,ordered):
    if(len(binaryvalue)==1):
        minterm=""
        if (binaryvalue[0] == 0):
            minterm += ordered[0]
        elif (binaryvalue[0] == 1):
            minterm += '(~' + ordered[0] + ')'
    else:
        minterm="("
        parantheses = 1
        for i in range(len(binaryvalue)):
            if(binaryvalue[i]==0):
                minterm+=ordered[i]
            elif(binaryvalue[i]==1):
                minterm+='(~'+ordered[i]+')'
            if(i<len(binaryvalue)-2):
                minterm+='&('
                parantheses+=1
            elif(i<len(binaryvalue)-1):
                minterm+='&'
        for j in range(parantheses):
            minterm+=')'
    return minterm
#end

vorodi=input()

#split (A>B) to A and B
c=0
index=0
for i in range(1,len(vorodi)):
    if(vorodi[i]=='('):
        c+=1
    elif(vorodi[i]==')'):
        c-=1
    if(c==0):
        index=i
        break
a=vorodi[1:index+1]
b=vorodi[index+2:len(vorodi)-1]

#extract atoms
atom_a=set()
atom_b=set()
atom_ab=set()
alaem={'(',')','|','>','~','&'}
for j in range(0,len(a)):
    if(a[j] not in alaem):
        atom_a.add(a[j])
for k in range(0,len(b)):
    if(b[k] not in alaem):
        atom_b.add(b[k])
atom_ab=atom_a.intersection(atom_b)

#ordered is a list of all of atoms that the same atoms come first of list
ordered=list()
for k in atom_ab:
    ordered.append(k)
allatom=(atom_a.union(atom_b))-atom_ab
for k in allatom:
    ordered.append(k)

#make the truth table of A
binarylist=[0 for k in range(len(ordered))]
truthtable_a=""
for i in range(2**(len(ordered))):
    if(answer(a,binarylist,ordered)):
        s="T"
    else:
        s="F"
    truthtable_a+=s
    generate(binarylist)

#make the truth table of A
binarylist=[0 for k in range(len(ordered))]
truthtable_b=""
for i in range(2**(len(ordered))):
    if(answer(b,binarylist,ordered)):
        s="T"
    else:
        s="F"
    truthtable_b+=s
    generate(binarylist)

#make the truth table of all of the possible interpolators
interpolator=generate_inter(truthtable_a,truthtable_b)

#this part check and find thr interpolators that can be pruduct with same atoms in A and B
finall_inter=set()
for some in interpolator:
    flag = True
    for j in range(0,len(some),2**len(allatom)):
        for k in range(j+1,j+(2**len(allatom))):
            if(some[j]!=some[k]):
                flag=False
                break
        if(not flag):
            break
    if(flag):
        finall_inter.add(some)

#this part reduce the truth table to pruduce with same atoms (example: TTTTFFTT -> TTFT)
finall_values=set()
for inter in finall_inter:
    tmp=""
    for j in range(0, len(inter), 2 ** len(allatom)):
        tmp+=inter[j]
    finall_values.add(tmp)

#this part convert truth table and pruduce formulas
inters=set()
for inter in finall_values:
    binaryvalue = [0 for i in range(len(atom_ab))]
    tmp=list()
    for x in inter:
        if(x=='T'):
            tmp.append(binary_to_minterm(binaryvalue,ordered))
        generate(binaryvalue)
    if(len(tmp)==0):
        inters.add("("+ordered[0]+"&(~"+ordered[0]+"))")
    elif(len(tmp)==1):
        inters.add(tmp[0])
    else:
        tmp2=""
        if(len(tmp)>1):
            tmp2+='('
        for i in range(len(tmp)-2):
            tmp2+=tmp[i]+"|("
        tmp2+=tmp[len(tmp)-2]+"|"+tmp[len(tmp)-1]
        for i in range(len(tmp)-1):
            tmp2+=')'
        inters.add(tmp2)

#print
for endgame in inters:
    print(endgame)