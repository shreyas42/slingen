'''
Created on Apr 18, 2012

@author: danieles
'''

import sys

from src.physical import Array, Scalars, Constant
#from math import *
from copy import copy

#from src.alexpr import *
from src.dsls.ll import Quantity, Operator, G, S, llLoop, llBlock, Sacc,\
    ParamMat, Assign
from src.dsls.sigmall import NewSum, Sum, llIf
from src.dsls.processing import computeIndependentSubexpr

  
class BindingTable(object):
    
    class myDict(dict):
        def __getitem__(self, *args, **kwargs):
            return dict.__getitem__(self, *args, **kwargs)
    
    def __init__(self):
        self.table = {}
#         self.table = BindingTable.myDict()
    
    def addBinding(self, matrix, physLayout):
        if matrix.name in self.table: return False
        self.table[matrix.name] = physLayout
        return True
    
    def add_binding_overwrite(self, matrix):
        if matrix.name in self.table: return False
        if matrix.attr['ow'] not in self.table:
            sys.exit("Could bind matrix %: %'s PhysLayout is missing." % (matrix.name, matrix.attr['ow']) )
        self.table[matrix.name] = self.table[matrix.attr['ow']]
        return True
        
    def delBindings(self, matrix):
        if matrix.name in self.table: del self.table[matrix.name] 
        
    def replacePhysicalLayout(self, oldPhys, newPhys):
        oldKeys = [ k for k, v in self.table.iteritems() if v == oldPhys]
        newPairs = [(k, newPhys) for k in oldKeys]
        self.table.update(newPairs)

    def replaceConnectedPhysicalLayout(self, oldPhys, newPhys, expr, stopOnNone=False):
        logic = expr.getOut().name
        phys = self.table[logic]
        replaced = set()
        if phys is None and not stopOnNone:
            if isinstance(expr, Operator):
                for sub in expr.inexpr:
                    #new change
                    stopOnNone = stopOnNone or isinstance(sub, ParamMat)
                    subPhys = self.table[sub.getOut().name]
                    replaced.update( self.replaceConnectedPhysicalLayout(subPhys, newPhys, sub, stopOnNone) )
        elif phys == oldPhys:
            self.table.update([(logic,newPhys)])
            replaced.update([phys])
            if isinstance(expr, Operator):
                for sub in expr.inexpr:
                    #new change
                    stopOnNone = stopOnNone or isinstance(sub, ParamMat)
                    replaced.update( self.replaceConnectedPhysicalLayout(oldPhys, newPhys, sub, stopOnNone) )
        return replaced
    
    def isBound(self, mat):
        return mat.name in self.table and self.table[mat.name] is not None # Any reason why the last condition should go??

    def existPhysicalLayout(self, phys):
        return phys in self.table.values()
            
    def getPhysicalLayout(self, matrix):
        return self.table[matrix.name]
    
    def resetTable(self):
        self.table.clear()
                
    def __str__(self):
        return str(sorted(self.table.items(), key=lambda x: x[0])) # str(self.table)

#bindingTable = BindingTable()

###################################################################################

class Binder(object):
    def __init__(self, context):
        super(Binder, self).__init__()
        self.context = context
    
    def apply(self, sllprog, opts):
        self._populate_binding_table(sllprog.mDict, opts) #mDict is the matrix information
        self._apply(sllprog.stmtList, opts) #apply the binding

    def _populate_binding_table(self, mat_dict, opts):
        import networkx as nx
        g = nx.DiGraph()
        for name,mat in mat_dict.iteritems():
            ow = mat.attr.get('ow', None)
            if ow is None:
                g.add_node(name)
            else:
                g.add_edge(name, ow)
        order_of_decl = nx.topological_sort(g, reverse=True)
        for name in order_of_decl:
            mat = mat_dict[name]
            if mat.attr.get('ow', None) is not None:
                self.context.bindingTable.add_binding_overwrite(mat)
            else:
                if not mat.attr['o'] and mat.isScalar():
                    physLayout = Scalars(mat.name, mat.size, opts, isIn=mat.attr['i'], isParam=True)
                else:
                    physLayout = Array(mat.name, mat.size, opts, isIn=mat.attr['i'], isOut=mat.attr['o'])
                if self.context.bindingTable.addBinding(mat, physLayout):
                    if mat.attr['t']:
                        physLayout.safelyScalarize = opts['scarep']
                        self.context.declare += [physLayout]
                    else:
                        self.context.signature += [physLayout]
        
        
    def _apply(self, expr, opts):
        if isinstance(expr, llBlock):
            for s in expr:
                self._apply(s, opts)
        elif isinstance(expr, llLoop):
            self._apply(expr.body, opts)
        elif isinstance(expr, llIf):
            for b in expr.bodys:
                self._apply(b, opts)
        else:
            getattr(self, expr.eq.__class__.__name__)(expr.eq, opts)

    def replaceConnectedPhysicalLayout(self, newPhys, expr, i):
        subPhys = self.context.bindingTable.getPhysicalLayout(expr.getInexprMat(i))
        self.context.bindingTable.replaceConnectedPhysicalLayout(subPhys, newPhys, expr.inexpr[i])
        if not self.context.bindingTable.existPhysicalLayout(subPhys):
            self.context.declare.remove(subPhys)
        
    def Assign(self, expr, opts):

        out = expr.getOut()
        if self.context.bindingTable.isBound(out):
            return
        
        getattr(self, expr.inexpr[0].__class__.__name__)(expr.inexpr[0], opts)
        getattr(self, expr.inexpr[1].__class__.__name__)(expr.inexpr[1], opts)

        mLhs = expr.getInexprMat(0)
        mRhs = expr.getInexprMat(1)
        
        if mLhs.attr['o']:
            src = mRhs
            dst = mLhs
            subExpr = expr.inexpr[1]
        else:
            src = mLhs
            dst = mRhs
            subExpr = expr.inexpr[0]
            
        #Replace the PhysLayout of the destination with the one of the source
#         if not subExpr.reqAss and not dst.attr['i']:
        if not subExpr.reqAss:
            srcPhys = self.context.bindingTable.getPhysicalLayout(src)
            dstPhys = self.context.bindingTable.getPhysicalLayout(dst)
            replaced = self.context.bindingTable.replaceConnectedPhysicalLayout(srcPhys, dstPhys, subExpr)
            for phys in replaced:
                if not self.context.bindingTable.existPhysicalLayout(phys):
                    self.context.declare.remove(phys)

    def Scalar(self, expr, opts):
        self.Matrix(expr, opts)
    
    def SquaredMatrix(self, expr, opts):
        self.Matrix(expr, opts)

    def Matrix(self, expr, opts):
        return
#         if self.context.bindingTable.isBound(expr):
#             return
#         if expr.attr.get('ow', None) is not None:
#             self.context.bindingTable.add_binding_overwrite(expr)
#         if not expr.attr['o'] and expr.isScalar():
#             physLayout = Scalars(expr.name, expr.size, opts, isIn=expr.attr['i'], isParam=True)
#         else:
#             physLayout = Array(expr.name, expr.size, opts, isIn=expr.attr['i'], isOut=expr.attr['o'])
#         if self.context.bindingTable.addBinding(expr, physLayout):
#             if expr.attr['t']:
#                 physLayout.safelyScalarize = opts['scarep']
#                 self.context.declare += [physLayout]
#             else:
#                 self.context.signature += [physLayout]
    
    def bindSimpleOp(self, expr, opts):
        out = expr.getOut()
        if self.context.bindingTable.isBound(out):
            return
        
        getattr(self, expr.inexpr[0].__class__.__name__)(expr.inexpr[0], opts)
        getattr(self, expr.inexpr[1].__class__.__name__)(expr.inexpr[1], opts)
        
        outPhys = Array(out.name, out.size, opts, safelyScalarize=opts['scarep'])
        if self.context.bindingTable.addBinding(out, outPhys):
            self.context.declare += [outPhys]
        
    def Add(self, expr, opts):
        self.bindSimpleOp(expr, opts)
        
    def Kro(self, expr, opts): # Temporarily only dealing with sca-mat mul
        self.bindSimpleOp(expr, opts)
            
    def Mul(self, expr, opts):
        self.bindSimpleOp(expr, opts)
    
    def PMul(self, expr, opts):
        self.bindSimpleOp(expr, opts)

    def bindSimpleUnary(self, expr, opts):
        out = expr.getOut()
        if self.context.bindingTable.isBound(out):
            return

        getattr(self, expr.inexpr[0].__class__.__name__)(expr.inexpr[0], opts)
        
        outPhys = Array(out.name, out.size, opts, safelyScalarize=opts['scarep'])
        if self.context.bindingTable.addBinding(out, outPhys):
            self.context.declare += [outPhys]

    def T(self, expr, opts):
        self.bindSimpleUnary(expr, opts)
        
    def HRed(self, expr, opts):
        out = expr.getOut()
        if self.context.bindingTable.isBound(out):
            return

        getattr(self, expr.inexpr[0].__class__.__name__)(expr.inexpr[0], opts)
        
        outPhys = Array(out.name, out.size, opts, safelyScalarize=opts['scarep'])
        if self.context.bindingTable.addBinding(out, outPhys):
            self.context.declare += [outPhys]
            
    def G(self, expr, opts):
        out = expr.getOut()
        if self.context.bindingTable.isBound(out):
            return

        getattr(self, expr.inexpr[0].__class__.__name__)(expr.inexpr[0], opts)
        
        subPhys = self.context.bindingTable.getPhysicalLayout(expr.getInexprMat(0))
        self.context.bindingTable.addBinding(out, subPhys)

    def S(self, expr, opts):

        out = expr.getOut()
        if self.context.bindingTable.isBound(out):
            return

        getattr(self, expr.inexpr[0].__class__.__name__)(expr.inexpr[0], opts)

        outPhys = Array(out.name, out.size, opts, safelyScalarize=opts['scarep'])
        if self.context.bindingTable.addBinding(out, outPhys):
            self.context.declare += [outPhys]
        if not isinstance(expr.inexpr[0], G):
            # If we directly scatter a gather we should keep the phys. layout separated
            # Otherwise we can bind the subexpr's phys. layout to a larger one where 
            # the op. S is supposed to scatter its input.
            self.replaceConnectedPhysicalLayout(outPhys, expr, 0)
    
    def Sum(self, expr, opts):
        
        out = expr.getOut()
        if self.context.bindingTable.isBound(out) or expr.isBinding:
            return
        
        expr.isBinding = True

        getattr(self, expr.inexpr[0].__class__.__name__)(expr.inexpr[0], opts)
        
        outPhys = None
        if isinstance(expr.inexpr[0], S):
            # In case of summing up scattered matrices we can either use S's physLayout
            # or making S use a new one. Here we go for the first option.
            outPhys = self.context.bindingTable.getPhysicalLayout(expr.getInexprMat(0))
            self.context.bindingTable.addBinding(out, outPhys)
        elif not isinstance(expr.inexpr[0], G):
            # If not directly summing up gathered matrices we can bind sub-expressions to the same
            # phys. layout.  
            outPhys = Array(out.name, out.size, opts, safelyScalarize=opts['scarep'])
            if self.context.bindingTable.addBinding(out, outPhys):
                self.context.declare += [outPhys]
            self.replaceConnectedPhysicalLayout(outPhys, expr, 0)
        
        # If the sum is accumulating, it should do it referring to same phys. layout
        if len(expr.inexpr) > 1:
            for i in range(1, len(expr.inexpr)):
                getattr(self, expr.inexpr[i].__class__.__name__)(expr.inexpr[i], opts)
                if not isinstance(expr.inexpr[i], G): 
                    # I'm not sure it can happen to have G subexprs when accumulating.
                    # By now checking such a condition just in case.
                    self.replaceConnectedPhysicalLayout(outPhys, expr, i)
        
        expr.isBinding = False

def old_detectScalarizableEquation(expr, context, explored, opts, prevSigma):
    '''
    When ScaRep is set many intermediate arrays can be completely removed (safelyScalarize).
    However, when this intermediate arrays are used to pass values between loop nests 
    they must be preserved.
    '''
    if isinstance(expr, Quantity): return

    ss = True # SafelyScalarize

    if (isinstance(expr, Sum) and not expr.isFullyUnrolled()) \
        or (expr in prevSigma and (len(prevSigma[expr]) != len(expr.pred) or not all(prevSigma[expr]))):
            ss = False

    out = expr.getOut()
    if context.bindingTable.isBound(out):
        phys = context.bindingTable.getPhysicalLayout(out)
        if isinstance(phys, Array):
            if not phys in context.signature: # Could be dropped - args are set not scalar. at binding time
#                 ss = False
                phys.safelyScalarize = ss
        
    if not expr in explored:
        explored.append(expr)
        subexprs = []
        isSum = False
        if isinstance(expr, Sum):
            isSum = True
            subexprs = computeIndependentSubexpr(expr.inexpr[0], expr.iList, explored, opts)
        else:
            subexprs = expr.inexpr
        
        for sub in subexprs:
            if isSum:
                if not sub in prevSigma:
                    prevSigma[sub] = [ ss ]
                else:
                    prevSigma[sub] += [ ss ]
            old_detectScalarizableEquation(sub, context, explored, opts, prevSigma)

def buildDUListEquation(expr, context, opts, duList, ctxList):

    if isinstance(expr, Quantity): return

    if isinstance(expr, NewSum):
        ctxList.append(expr)
        buildDUListEquation(expr.inexpr[0], context, opts, duList, ctxList)
        ctxList.pop()
    else:
        if isinstance(expr, Assign):
            lhs = expr.inexpr[0]
            if lhs.attr['t']:
                if context.bindingTable.isBound(lhs):
                    phys = context.bindingTable.getPhysicalLayout(lhs)
                    duList.append([phys])
        if isinstance(expr, ParamMat):
            out = expr.getOut()
            if context.bindingTable.isBound(out):
                phys = context.bindingTable.getPhysicalLayout(out)
                dus = filter(lambda du: du[0] == phys, duList)
                if dus:
                    ctxl = [ ctx for ctx in ctxList if ctx.idx in expr.fL.func or ctx.idx in expr.fR.func ]
                    if isinstance(expr, Sacc):
                        du,top = (dus[-2],dus[-1]) if len(dus[-1]) == 1 else (dus[-1], None)
                        du.append( ('u', ctxl) )
                        if top is None:
                            newdu = [ du[0], ('d', ctxl) ] 
                            duList.append(newdu)
                        else:
                            top.append( ('d', ctxl) )
                    elif isinstance(expr, G):
                        du = dus[-2] if len(dus[-1]) == 1 else dus[-1]
                        du.append( ('u', ctxl) )
                    elif isinstance(expr, S):
                        dus[-1].append( ('d', ctxl) )
        for sub in expr.inexpr:
            buildDUListEquation(sub, context, opts, duList, ctxList)

def buildDUList(expr, context, opts, duList, ctxList=None):
    ctxList = [] if ctxList is None else ctxList
    if isinstance(expr, llBlock):
        for s in expr:
            buildDUList(s, context, opts, duList, ctxList)
    elif isinstance(expr, llLoop):
        ctxList.append(expr)
        buildDUList(expr.body, context, opts, duList, ctxList)
        ctxList.pop()
    elif isinstance(expr, llIf):
        for b in expr.bodys:
            buildDUList(b, context, opts, duList, ctxList)
    else:
        buildDUListEquation(expr.eq, context, opts, duList, ctxList)
    
def detectNonScalarizableMats(duList, opts):
    nonScal = []
    for du in duList: # [ MATLayout, ('d', [for0, for1, ...]), ('u', [for0, for1, ...]), ... ]
        if du[0] in nonScal:
            continue
        dctxs = du[1][1]
        for u in du[2:]:
            i=0
            while i < min(len(dctxs),len(u[1])):
#                 if dctxs[i].idx != u[1][i].idx:
                if id(dctxs[i]) != id(u[1][i]):
                    break
                i += 1
            checkUnrollable = dctxs[i:] + u[1][i:]
            for f in checkUnrollable:
                L = f.ub-f.lb+1
                if (opts['unroll'][str(f.idx)] == 0) or ((opts['unroll'][str(f.idx)] > 0) and not L.is_Number):
                    du[0].safelyScalarize = False
                    nonScal.append(du[0])
                    break

def bindExpression(sllprog, context, opts=None):
#     binder = Binder(context)
#    binder = NewBinder(context)
    binder = EnhancedBinder(context)
    binder.apply(sllprog, opts)
    tsig = []
    for m in opts['inoutorder']:
        tsig += filter(lambda phys: phys.name == m, context.signature)
    del context.signature[:]
    context.signature = tsig
    if opts is not None and opts.get('scarep', False):
        duList = []
        buildDUList(sllprog.stmtList, context, opts, duList)
        detectNonScalarizableMats(duList, opts)


class NewBinder(Binder):
    def __init__(self, context):
        super(NewBinder, self).__init__(context)
        self.scattering = 0

    def G(self, expr, opts):
        out = expr.getOut()
        if out.size[0] <= opts['nu'] and out.size[1] <= opts['nu']: 
            getattr(self, expr.inexpr[0].__class__.__name__)(expr.inexpr[0], opts)
            
            subPhys = self.context.bindingTable.getPhysicalLayout(expr.getInexprMat(0))
            self.context.bindingTable.addBinding(out, subPhys)
            out.fL, out.fR = expr.fL, expr.fR
        else:
            if self.context.bindingTable.isBound(out):
                return
            getattr(self, expr.inexpr[0].__class__.__name__)(expr.inexpr[0], opts)
            outPhys = Array(out.name, out.size, opts)
            if self.context.bindingTable.addBinding(out, outPhys):
                self.context.declare += [outPhys]
            
    def S(self, expr, opts):
        self.scattering += 1
        out = expr.getOut()
        if self.context.bindingTable.isBound(out):
            return
        getattr(self, expr.inexpr[0].__class__.__name__)(expr.inexpr[0], opts)

        sub = expr.getInexprMat(0)
        safelyScalarize = opts['scarep'] and sub.size[0] <= opts['nu'] and sub.size[1] <= opts['nu']
        outPhys = Array(out.name, out.size, opts, safelyScalarize=safelyScalarize)
        if self.context.bindingTable.addBinding(out, outPhys):
            self.context.declare += [outPhys]
####### PhysLayout replacement should be delayed until structures and accesses are computed ############
 
#         if sub.size[0] <= opts['nu'] and sub.size[1] <= opts['nu']: 
#             outPhys = Array(out.name, out.size, opts, safelyScalarize=opts['scarep'])
#             if self.context.bindingTable.addBinding(out, outPhys):
#                 self.context.declare += [outPhys]
#             
#             if not isinstance(expr.inexpr[0], G):
#                 # If we directly scatter a gather we should keep the phys. layout separated
#                 # Otherwise we can bind the subexpr's phys. layout to a larger one where 
#                 # the op. S is supposed to scatter its input.
#                 self.replaceConnectedPhysicalLayout(outPhys, expr, 0)
#                 sub.fL, sub.fR = expr.fL, expr.fR
#         else:
#             outPhys = Array(out.name, out.size, opts)
#             if self.context.bindingTable.addBinding(out, outPhys):
#                 self.context.declare += [outPhys]
            
        self.scattering -= 1

    def Sacc(self, expr, opts):
#         self.S(expr, opts)
        self.scattering += 1
        out = expr.getOut()
        if self.context.bindingTable.isBound(out):
            return
        print(expr.inexpr[0].__class__.__name__)
        getattr(self, expr.inexpr[0].__class__.__name__)(expr.inexpr[0], opts)

        #added code
        #end of added code

        sub = expr.getInexprMat(0)

        if sub.size[0]*sub.size[1] <= opts['nu']*opts['nu']:
            outPhys = Array(out.name, out.size, opts, safelyScalarize=opts['scarep'])
            if self.context.bindingTable.addBinding(out, outPhys):
                self.context.declare += [outPhys]
        else:
            outPhys = Array(out.name, out.size, opts)
            if self.context.bindingTable.addBinding(out, outPhys):
                self.context.declare += [outPhys]

        self.scattering -= 1

    def Neg(self, expr, opts):
        self.bindSimpleUnary(expr, opts)

    def NewSum(self, expr, opts):
        out = expr.getOut()
        if self.context.bindingTable.isBound(out):
            return
        
        getattr(self, expr.inexpr[0].__class__.__name__)(expr.inexpr[0], opts)
        
        outPhys = self.context.bindingTable.getPhysicalLayout(expr.getInexprMat(0))
        self.context.bindingTable.addBinding(out, outPhys)

    #modifying the Add() and the Mul() methods

    def Add(self, expr, opts):
        out = expr.getOut()
        if self.context.bindingTable.isBound(out):
            return
        getattr(self, expr.inexpr[0].__class__.__name__)(expr.inexpr[0], opts)
        getattr(self, expr.inexpr[1].__class__.__name__)(expr.inexpr[1], opts)
        #added code
        #out1 = expr.inexpr[0].getOut()
        #out2 = expr.inexpr[1].getOut()
        #if not self.context.bindingTable.isBound(out1) or not self.context.bindingTable.isBound(out2):
        #    print('FMA detected')
        #end of added code
        if not self.scattering:
            self.context.bindingTable.addBinding(out, None)
        else:
            self.bindSimpleOp(expr, opts)

#added the overridden Mul() method

    def Mul(self , expr , opts):
        #here's what needs to happen
        #first off we gotta check if the parent of the expression is a + operator
        #print('Overriden Mul() method ')
        out  = expr.getOut()
        if self.context.bindingTable.isBound(out):
            return
        parent = expr.pred[0]
        #print(parent[0].__class__.__name__)
        if parent[0].__class__.__name__ in [ 'Add' , 'Sacc' ]:
        #    print('Some dummy print statement')
            self.context.bindingTable.addBinding(out , None)
            #okay so it can detect if its parent is a add expression
            #so what needs to happen next is that it should not be bound
            getattr(self , expr.inexpr[0].__class__.__name__)(expr.inexpr[0] , opts)
            getattr(self , expr.inexpr[1].__class__.__name__)(expr.inexpr[1] , opts)
        else:
            self.bindSimpleOp(expr , opts)
#end of overridden method

    def LDiv(self, expr, opts):
        self.bindSimpleOp(expr, opts)

    def Div(self, expr, opts):
        self.bindSimpleOp(expr, opts)

    def Sqrt(self, expr, opts):
        self.bindSimpleUnary(expr, opts)

    def Sub(self, expr, opts):
        self.bindSimpleOp(expr, opts)

    def Iv(self, expr, opts):        
        out = expr.getOut()
        if self.context.bindingTable.isBound(out):
            return

        getattr(self, expr.inexpr[0].__class__.__name__)(expr.inexpr[0], opts)
        
        outPhys = self.context.bindingTable.getPhysicalLayout(expr.getInexprMat(0))
        self.context.bindingTable.addBinding(out, outPhys)

    def LowerTriangular(self, expr, opts):
        self.Matrix(expr, opts)

    def LowerUnitTriangular(self, expr, opts):
        self.Matrix(expr, opts)

    def UpperTriangular(self, expr, opts):
        self.Matrix(expr, opts)

    def UpperUnitTriangular(self, expr, opts):
        self.Matrix(expr, opts)

    def Symmetric(self, expr, opts):
        self.Matrix(expr, opts)
    
    def ConstantMatrix(self, expr, opts):
        
        if self.context.bindingTable.isBound(expr):
            return
        
        physLayout = Constant()
        self.context.bindingTable.addBinding(expr, physLayout)

    def IdentityMatrix(self, expr, opts):
        self.ConstantMatrix(expr, opts)

    def AllEntriesConstantMatrixWithValue(self, expr, opts):
        self.ConstantMatrix(expr, opts)

#MARKER FOR CHANGE IN CODE
#assuming that we arent introducing a new physical layout and are using the same concept of an array
#the only change is what we pass as size to the Array constructor?
#do any functions(outside of the Binder classes) need to change?

class EnhancedBinder(NewBinder):
    def __init__(self , context):
        super(EnhancedBinder , self).__init__(context)

    def _populate_binding_table(self, mat_dict, opts):
        import networkx as nx
        g = nx.DiGraph()
        for name,mat in mat_dict.iteritems():
            ow = mat.attr.get('ow', None)
            if ow is None:
                g.add_node(name)
            else:
                g.add_edge(name, ow)
        order_of_decl = nx.topological_sort(g, reverse=True)
        for name in order_of_decl:
            mat = mat_dict[name]
            isComplex = False
            #phys_size = list(mat.size)
            if mat.get_field() == 'complex':
            #    phys_size[1] *= 2
                isComplex = True
            #phys_size = tuple(phys_size)

            if mat.attr.get('fieldinfo' , None) is not None and len(mat.attr['fieldinfo']) > 1:
                fieldParam = mat.attr['fieldinfo'][1]
            else:
                fieldParam = None

            if mat.attr.get('ow', None) is not None:
                self.context.bindingTable.add_binding_overwrite(mat)
            else:
                if not mat.attr['o'] and mat.isScalar():
                    physLayout = Scalars(mat.name, mat.size, opts, isIn=mat.attr['i'], isParam=True)
                else:
                    #physLayout = Array(mat.name, phys_size, opts, isIn=mat.attr['i'], isOut=mat.attr['o'] , field = fieldParam)
                    physLayout = Reference.createArray(mat , opts , useIn = True , useOut = True , useComplex = isComplex , field = fieldParam)
                if self.context.bindingTable.addBinding(mat, physLayout):
                    if mat.attr['t']:
                        physLayout.safelyScalarize = opts['scarep']
                        self.context.declare += [physLayout]
                    else:
                        self.context.signature += [physLayout]

    def bindSimpleOp(self , expr , opts):
        out = expr.getOut()
        if self.context.bindingTable.isBound(out):
            return

        getattr(self, expr.inexpr[0].__class__.__name__)(expr.inexpr[0], opts)
        getattr(self, expr.inexpr[1].__class__.__name__)(expr.inexpr[1], opts)
        is_complex = False
        #phys_size = list(out.size)
        if out.get_field() == 'complex':
            #phys_size[1] *= 2
            is_complex = True
        #phys_size = tuple(phys_size)
        #using the hackish way to get shit done , will have to change later
        outPhys = Reference.createArray(out , opts , useScalarize=True , useComplex = is_complex)
        '''if is_complex:
            outPhys = Array(out.name, phys_size, opts, safelyScalarize=opts['scarep'] , field = 'BlkInterLeaved')
        else:
            outPhys = Array(out.name, phys_size, opts, safelyScalarize=opts['scarep'])
        '''
        if self.context.bindingTable.addBinding(out, outPhys):
            self.context.declare += [outPhys]

    def bindSimpleUnary(self, expr, opts):
        out = expr.getOut()
        if self.context.bindingTable.isBound(out):
            return

        getattr(self, expr.inexpr[0].__class__.__name__)(expr.inexpr[0], opts)
        #phys_size = list(out.size)
        is_complex = False
        if out.get_field() == 'complex':
        #    phys_size[1] *= 2
            is_complex = True
        #phys_size = tuple(phys_size)

        outPhys = Reference.createArray(out , opts , useScalarize=True , useComplex = is_complex)
        '''if is_complex:
            outPhys = Array(out.name, phys_size, opts, safelyScalarize=opts['scarep'] , field = 'BlkInterLeaved')
        else:
            outPhys = Array(out.name, phys_size, opts, safelyScalarize=opts['scarep'])
        '''
        if self.context.bindingTable.addBinding(out, outPhys):
            self.context.declare += [outPhys]

    def G(self , expr , opts):
        out = expr.getOut()
        if out.size[0] <= opts['nu'] and out.size[1] <= opts['nu']:
            getattr(self, expr.inexpr[0].__class__.__name__)(expr.inexpr[0], opts)

            subPhys = self.context.bindingTable.getPhysicalLayout(expr.getInexprMat(0))
            self.context.bindingTable.addBinding(out, subPhys)
            out.fL, out.fR = expr.fL, expr.fR
        else:
            if self.context.bindingTable.isBound(out):
                return
            getattr(self, expr.inexpr[0].__class__.__name__)(expr.inexpr[0], opts)
            #phys_size = list(out.size)
            is_complex = False
            if out.get_field() == 'complex':
            #    phys_size[1] *= 2
                is_complex = True
            #phys_size = tuple(phys_size)
            outPhys = Reference.createArray(out , opts , useScalarize=True , useComplex = is_complex)
            '''
            if is_complex:
                outPhys = Array(out.name, phys_size, opts , field = 'BlkInterLeaved')
            else:
                outPhys = Array(out.name, phys_size, opts)
            '''
            if self.context.bindingTable.addBinding(out, outPhys):
                self.context.declare += [outPhys]

    def S(self, expr, opts):
        self.scattering += 1
        out = expr.getOut()
        if self.context.bindingTable.isBound(out):
            return
        getattr(self, expr.inexpr[0].__class__.__name__)(expr.inexpr[0], opts)

        sub = expr.getInexprMat(0)
        safelyScalarize = opts['scarep'] and sub.size[0] <= opts['nu'] and sub.size[1] <= opts['nu']
        #phys_size = list(out.size)
        is_complex = False
        if out.get_field() == 'complex':
        #    phys_size[1] *= 2
            is_complex = True
        #phys_size = tuple(phys_size)
        outPhys = Reference.createArray(out , opts , useScalarize=True , useComplex = is_complex)
        '''
        if is_complex:
            outPhys = Array(out.name, phys_size, opts, safelyScalarize=safelyScalarize , field = 'BlkInterLeaved')
        else:
            outPhys = Array(out.name, phys_size, opts, safelyScalarize=safelyScalarize)
        '''
        if self.context.bindingTable.addBinding(out, outPhys):
            self.context.declare += [outPhys]

        self.scattering -= 1

    def Sacc(self, expr, opts):
        #         self.S(expr, opts)
        self.scattering += 1
        out = expr.getOut()
        if self.context.bindingTable.isBound(out):
            return
        print(expr.inexpr[0].__class__.__name__)
        getattr(self, expr.inexpr[0].__class__.__name__)(expr.inexpr[0], opts)

        sub = expr.getInexprMat(0)

        #phys_size = list(out.size)
        is_complex = False
        if out.get_field() == 'complex':
      #      phys_size[1] *= 2
            is_complex = True
       # phys_size = tuple(phys_size)

        if sub.size[0]*sub.size[1] <= opts['nu']*opts['nu']:

            outPhys = Reference.createArray(out , opts , useScalarize=True , useComplex = is_complex)
            '''if is_complex:
                outPhys = Array(out.name, phys_size, opts, safelyScalarize=opts['scarep'] , field = 'BlkInterLeaved')
            else:
                outPhys = Array(out.name, phys_size, opts, safelyScalarize=opts['scarep'])
            '''
            if self.context.bindingTable.addBinding(out, outPhys):
                self.context.declare += [outPhys]
        else:

            outPhys = Reference.createArray(out , opts , useComplex = is_complex)
            '''if is_complex:
                outPhys = Array(out.name, phys_size, opts , field = 'BlkInterLeaved')
            else:
                outPhys = Array(out.name, phys_size, opts)
            '''
            if self.context.bindingTable.addBinding(out, outPhys):
                self.context.declare += [outPhys]

        self.scattering -= 1

###################################################################################

class Reference(object):
    def __init__(self, matrix, physLayout):
        self.matrix = matrix
        self.physLayout = physLayout

    @staticmethod
    def createArray(mat , opts , useIn = False ,useOut = False , useScalarize =False, useComplex=False , field = None):
        if useComplex:
            phys_size = list(mat.size)
            if opts['complexlayout'] == 'BlkInterLeaved':
                #assuming that nu is going to be the block size
                blocksize = opts['nu']
                if phys_size[1] % blocksize != 0:
                    phys_size[1] += (blocksize - (phys_size[1] % blocksize))
                phys_size[1] *= 2

            elif opts['complexlayout'] == 'Split':
                pass

            phys_size = tuple(phys_size)

            if useIn and useOut:
                if field is None:
                    return Array(mat.name, phys_size, opts, isIn=mat.attr['i'] ,isOut = mat.attr['o'], field = opts['complexlayout'] )
                else:
                    return Array(mat.name, phys_size, opts, isIn=mat.attr['i'] ,isOut = mat.attr['o'], field = field )

            elif useScalarize:
                if field is None:
                    return Array(mat.name, phys_size, opts, safelyScalarize=opts['scarep'], field = opts['complexlayout'] )
                else:
                    return Array(mat.name, phys_size, opts, safelyScalarize=opts['scarep'], field = field )

            else:
                if field is None:
                    return Array(mat.name, phys_size, opts, field = opts['complexlayout'] )
                else:
                    return Array(mat.name, phys_size, opts, field = field )

        else:
            if useIn and useOut:
                return Array(mat.name, mat.size, opts, isIn=mat.attr['i'] ,isOut = mat.attr['o'] )
            elif useScalarize:
                return Array(mat.name, mat.size, opts, safelyScalarize=opts['scarep'] )
            else :
                return Array(mat.name, mat.size, opts)

    @staticmethod
    def whatRef(PhysLayout):
        if isinstance(PhysLayout, Array):
            if PhysLayout.getField() == 'BlkInterLeaved':
                return RowMajorInterArrayReference
            elif PhysLayout.getField() == 'Split':
                pass
            else:
                return ArrayReference
        if isinstance(PhysLayout, Scalars):
            return ScalarsReference
        if isinstance(PhysLayout, Constant):
            return ConstantReference
        else:
            return None

    def __eq__(self, other):
        return self.physLayout == other.physLayout

    def __str__(self):
        return self.matrix.name + " -> " + self.physLayout.name
    
class ConstantReference(Reference):
    def __init__(self, matrix, physLayout):
        super(ConstantReference, self).__init__(matrix, physLayout)

class ExplicitPhysicalReference(Reference):
    def __init__(self, matrix, physLayout):
        super(ExplicitPhysicalReference, self).__init__(matrix, physLayout)
    
    def getLinIdx(self, key):
        return None

    def __getitem__(self, key):
        linIdx = self.getLinIdx(key)
        return self.physLayout[linIdx]

#associates matrices(mathematical concept) to physical storage
class ArrayReference(ExplicitPhysicalReference):
    def __init__(self, matrix, physLayout):
        super(ArrayReference, self).__init__(matrix, physLayout)


    #mathematically we access element at i,j : this function maps this to a linear index in the physical storage of this matrix
    def getLinIdx(self, key):
        idx = copy(self.matrix.getOrigin()) #get the index of the origin - again this is a logical concept
        #when would the origin not be 0,0?
        #cause matrix is just a logical concept - so it should be that the origin is at 0,0
        idx[0] += key[0]
        idx[1] += key[1]
        
        return idx[0]*self.physLayout.pitch + idx[1] #in any event this returns the linear index that the element is at in the physical storage(array)
        
    def pointerAt(self, key):
#         idx = copy(self.matrix.getOrigin())
#         idx[0] += key[0]
#         idx[1] += key[1]
#         
#         linIdx = idx[0]*self.physLayout.pitch + idx[1]
        linIdx = self.getLinIdx(key)
        return self.physLayout.pointerAt(linIdx)

    def isCorner(self, key):
        linIdx = self.getLinIdx(key)
        return self.physLayout.size-1 == linIdx
#     def __getitem__(self, key):
#         idx = copy(self.matrix.getOrigin())
#         idx[0] += key[0]
#         idx[1] += key[1]
#         
#         linIdx = idx[0]*self.physLayout.pitch + idx[1]
#         linIdx = self.getLinIdx(key)
#         return self.physLayout[linIdx]

#design choice : should this extend ArrayReference or ExplicitPhysicalReference
class RowMajorInterArrayReference(ExplicitPhysicalReference):
    def __init__(self , matrix , physLayout):
        super(RowMajorInterArrayReference , self).__init__(matrix , physLayout)

    #the new getLinIdx function
    def getLinIdx(self , key):
        import sympy
        idx = copy(self.matrix.getOrigin())
        idx[0] += key[0]
        idx[1] += key[1]
        selection = key[2] # 0 for the real part , 1 for the imaginary part
        block_size = self.physLayout.blocksize #this is under consideration
        #this statement is problematic
        block_num = sympy.floor(idx[1] / block_size)
        e = idx[1] % block_size

        #if self.matrix.size[1] % block_size == 0:
        #    x = block_size
        #else:
        #    x = self.matrix.size[1] % block_size
        return (idx[0] * self.physLayout.pitch) + (2 * block_size * block_num) + e + (selection * block_size)
        #return (idx[0] * self.physLayout.pitch) + (2 * block_size * block_num) + e + (selection * x)

    def pointerAt(self , key):
        linIdx = self.getLinIdx(key)
        return self.physLayout.pointerAt(linIdx)

    def isCorner(self , key):
        linIdx = self.getLinIdx(key)
        return self.physLayout.size-1 == linIdx

class ScalarsReference(ExplicitPhysicalReference):
    def __init__(self, matrix, physLayout):
        super(ScalarsReference, self).__init__(matrix, physLayout)

    def getLinIdx(self, key):
        #this method will need to change at some point
        idx = copy(self.matrix.getOrigin())
        idx[0] += key[0]
        idx[1] += key[1]
        
        return idx[0]*self.physLayout.size[1] + idx[1]

    def isCorner(self, key):
        return True 

#     def __getitem__(self, key):
#         idx = copy(self.matrix.getOrigin())
#         idx[0] += key[0]
#         idx[1] += key[1]
#         
#         linIdx = idx[0]*self.physLayout.size[1] + idx[1]
#         linIdx = self.getLinIdx(key)
#         return self.physLayout[linIdx]


###################################################################################


def getReference(context, matrix):
    '''Get physical layout reference.'''
    if not context.bindingTable.isBound(matrix): #is the matrix bound to some physical layout?
        return None
    physLayout = context.bindingTable.getPhysicalLayout(matrix) # retrieve its physical layout #if i need to add complex matrices then the matrix needs to be of a different physical layout
    #changed from a static method - the question is why exactly? Didn't really get that
    #you could still make it static and pass the matrix as a parameter
    #this way you could query the field of the matrix and return the right Reference type
    Ref = Reference.whatRef(physLayout) #what kind of reference object will it be?
    return Ref(matrix, physLayout) #return an object of the correct reference type passing arguments: matrix and the retrieved physical layout

###################################################################################


if __name__ == "__main__":
    pass    
