import sys,os,re,time


globalvars = []
hfiles = "after_opts.h all_atom.h misc.h param.h fullatom_energies.h"
globalvarfinder = re.compile("extern\s+(\S+)\s+(\S+?);",re.M)
for hfile in hfiles.split():
    h = open(hfile).read()
    g = globalvarfinder.findall(h)
    for v in g:
        globalvars.append((v[1],v[0],hfile))
    
#print globalvars

def pad(strings,left=True):
    l = max([len(s) for s in strings])
    if left:
        return [s+" "*(l-len(s)) for s in strings]
    else:
        return [" "*(l-len(s))+s for s in strings]
    
def makecomment(rettype, ns, fname, args, body, other={}):
    commenttemplate = """
////////////////////////////////////////////////////////////////////// 
/// @begin %(fname)s
///
/// @brief 
%(brief)s
/// 
/// @detailed 
%(detailed)s
/// 
%(params)s
/// 
/// @return
%(ret)s
/// 
/// @global_read
%(globalread)s
/// 
/// @global_write
%(globalwrite)s
/// 
/// @remarks 
%(remarks)s
/// 
/// @references 
%(references)s
/// 
/// @authors %(author)s
/// 
/// @last_modified %(date)s
/////////////////////////////////////////////////////////////////////////
"""

    remarks = '/// TODO: fill me in!'    
    if ns != "":
        if ns[-2:] == "::":
            ns = ns[:-2]
        remarks = "/// TODO: fill me in!\n/// In namespace %s"%ns

    if args.strip() == "":
        params = "/// No Params"
    else:
        args = args.split(',')
        def kind(p):
            if ("&" in p or "*" in p) and not "const" in p: return '[in,out?]'
            else:                                           return "[in]"
        ptypes = pad([" ".join(p.split()[:-1]) for p in args],left=True)
        pnames = pad([         p.split()[ -1]  for p in args],left=True)
        pkinds = pad([    kind(p)              for p in args],left=True)
        params = ""
        for ii in range(len(args)):
            p     = args[ii]
            ptype = ptypes[ii]
            pname = pnames[ii]
            pkind = pkinds[ii]
            params += "/// @param%s %s %s - TODO: fill me in!\n"%(pkind,ptype,pname)
        params = params[:-1]

    body = re.sub("//.*;","",body)
    body = re.sub("/\*.*\*/;","",body,re.M)
    rg = [s for s in globalvars if s[0] in body and containsVar(var=s[0],body=body)]
    wg = [s for s in globalvars if s[0] in body and   writesVar(var=s[0],body=body)]
    #rg = [("varname:","var type:","header:")] + rg
    #wg = [("varname:","var type:","header:")] + wg
    globalread = '/// TODO: fill me in!'
    if len(rg) > 1:
        globalread += ", but here are some guesses:"
        grname = pad([r[0] for r in rg])
        grtype = pad([r[1] for r in rg])
        grfile = pad([r[2] for r in rg])
        for ii in range(len(rg)):
            globalread += "\n///     %s %s in %s"%(grtype[ii],grname[ii],grfile[ii])
    globalwrite = '/// TODO: fill me in!'
    if len(wg) > 1:
        globalwrite += ", but here are some guesses:"
        grname = pad([r[0] for r in wg])
        grtype = pad([r[1] for r in wg],left=False)
        grfile = pad([r[2] for r in wg])
        for ii in range(len(wg)):
            globalwrite += "\n///     %s %s in %s"%(grtype[ii],grname[ii],grfile[ii])
        
    fields = {
        'fname':fname,
        'brief':'/// TODO: fill me in!',
        'detailed':'/// TODO: fill me in!',
        'params':params,
        'ret':'/// TODO: fill me in!\n/// '+rettype,
        'globalread':globalread,
        'globalwrite':globalwrite,
        'remarks':remarks,
        'references':'/// TODO: fill me in!',
        'author': os.popen2("whoami")[1].read()[:-1],
        'date': time.ctime()
        }
    for k in other:
        fields[k] = other[k]
    return commenttemplate%fields

def containsVar(var,body):
    begin = "[ \t\n\r\f\v,;(){}=+\-/*&^!%|?:\[\]]"
    end   = "[ \t\n\r\f\v,;(){}=+\-/*&^!%|?:\[\].]"
    tmp = re.compile(begin+var+end,re.M).search(body)
    if tmp is not None: return True
    return False

def writesVar(var,body):
    begin = "[ \t\n\r\f\v,;(){}=+\-/*&^!%|?:\[\]]"
    end   = "\s*(\(.*?\))?(\[.*?\])?\s*[+-/*]?=[^=]"
    tmp = re.compile(begin+var+end,re.M).search(body)
    if tmp is not None: return True
    begin = "[ \t\n\r\f\v,;(){}=+\-/*&^!%|?:\[\]]"
    end   = "[ \t\n\r\f\v,;(){}=+\-/*&^!%|?:\[\].]"
    tmp = re.compile(begin    +var+"\s*\+\+",re.M).search(body)
    if tmp is not None: return True
    tmp = re.compile("\+\+\s*"+var+end      ,re.M).search(body)
    if tmp is not None: return True
    return False
    
retype1 = "^([^/#;{}() \t\n\r\f\v]+?)\s+"
retype2 = "^"
renmsp1 = "((?:\S*?::)*)"
renmsp2 = "((?:\S*?::)+)"
rename = "([^/#;{}() \t\n\r\f\v]+?)\s*?"
reargs = "\((.*?)\).*?"
rebody = "({[^{]*}|{.*?^})"
regexp1 = retype1+renmsp1+rename+reargs+rebody
regexp2 = retype2+renmsp2+rename+reargs+rebody
regexp  = regexp1+'|'+regexp2

f = open(sys.argv[1])
code  = f.read()
f.close()

prev = 0
for match in re.compile(regexp,re.DOTALL|re.M).finditer(code):
    start = match.start()
    end   = match.end()
    span  = match.span()
    if match.group(1) is not None:
        ret  = match.group(1)
        ns   = match.group(2)
        name = match.group(3)
        args = match.group(4)
        body = match.group(5)
    else:
        ret  = "None (constructor)"
        ns   = match.group(6)
        name = match.group(7)
        args = match.group(8)
        body = match.group(9)
    stub = makecomment(ret,ns,name,args,body)
    print code[prev:start] 
    print stub
    prev = start
print code[prev:]
    
