import sys
from sympy import sympify

from islpy import Set, Map

from src.dsls.ll import Matrix, ZeroMatrix, Symmetric, LowerTriangular, UpperTriangular, LowerUnitTriangular, UpperUnitTriangular, IdentityMatrix, AllEntriesConstantMatrix

from src.binding import getReference, ScalarsReference
from src.irbase import RValue, Pointer, VecAccess, VecDest, MovStatement, Mov, Comment, AddressOf, sa, V, DebugPrint, icode
from src.isas.isabase import ISA, Loader, Storer, LoadReplacer

#which instructions are we going to include
#considering only double-precision floating point numbers for now

class mm256FmaddPd(RValue):
    def __init__(self,src0,src1,src2):
        super(mm256FmaddPd,self).__init__()
        self.srcs += [ src0 , src1 , src2 ]

        def unparse(self,indent):
            return indent + "_mm256_fmadd_pd(" + self.srcs[0].unparse("") + ", " + self.srcs[1].unparse("") + ", " + self.srcs[2].unparse("") + " )"

        def printInst(self,indent):
            return indent + "mm256FmaddPd( " + self.srcs[0].unparse("") + ", " + self.srcs[1].unparse("") + ", " + self.srcs[2].unparse("") + " )"

class mm256FmaddsubPd(Rvalue):
    def __init__(self,src0,src1,src2):
        super(mm256FmaddsubPd,self).__init__()
        self.srcs += [ src0 , src1 , src2 ]

    def unparse(self,indent):
        return indent + "_mm256_fmaddsub_pd(" + self.srcs[0].unparse("") + ", " + self.srcs[1].unparse("") + ", " + self.srcs[2].unparse("") + " )"

    def printInst(self,indent):
        return indent + "mm256FmaddsubPd( " + self.srcs[0].unparse("") + ", " + self.srcs[1].unparse("") + ", " + self.srcs[2].unparse("") + " )"

class mm256FmsubPd(Rvalue):
    def __init__(self,src0,src1,src2):
        super(mm256FmsubPd,self).__init__()
        self.srcs += [ src0 , src1 , src2 ]

    def unparse(self,indent):
        return indent + "_mm256_fmsub_pd(" + self.srcs[0].unparse("") + ", " + self.srcs[1].unparse("") + ", " + self.srcs[2].unparse("") + " )"

    def printInst(self,indent):
        return indent + "mm256FmsubPd( " + self.srcs[0].unparse("") + ", " + self.srcs[1].unparse("") + ", " + self.srcs[2].unparse("") + " )"

class mm256FmsubaddPd(Rvalue):
    def __init__(self,src0,src1,src2):
        super(mm256FmsubaddPd,self).__init__()
        self.srcs += [ src0 , src1 , src2 ]

    def unparse(self,indent):
        return indent + "_mm256_fmsubadd_pd(" + self.srcs[0].unparse("") + ", " + self.srcs[1].unparse("") + ", " + self.srcs[2].unparse("") + " )"

    def printInst(self,indent):
        return indent + "mm256FmsubaddPd( " + self.srcs[0].unparse("") + ", " + self.srcs[1].unparse("") + ", " + self.srcs[2].unparse("") + " )"

class mm256FnmaddPd(Rvalue):
    def __init__(self,src0,src1,src2):
        super(mm256FnmaddPd,self).__init__()
        self.srcs += [ src0 , src1 , src2 ]

    def unparse(self,indent):
        return indent + "_mm256_fnmadd_pd(" + self.srcs[0].unparse("") + ", " + self.srcs[1].unparse("") + self.srcs[2].unparse("") + " )"

    def printInst(self,indent):
        return indent + "mm256FnmaddPd(" + self.srcs[0].unparse("") + ", " + self.srcs[1].unparse("") + self.srcs[2].unparse("") + " )"

class mm256FnmsubPd(Rvalue):
    def __init__(self,src0,src1,src2):
        super(mm256FnmsubPd,self).__init__()
        self.srcs += [ src0 , src1 , src2 ]

    def unparse(self,indent):
        return indent + "_mm256_fnmsub_pd(" + self.srcs[0].unparse("") + ", " + self.srcs[1].unparse("") + self.srcs[2].unparse("") + " )"

    def printInst(self,indent):
        return indent + "mm256FnmsubPd(" + self.srcs[0].unparse("") + ", " + self.srcs[1].unparse("") + self.srcs[2].unparse("") + " )"

class _Fma4BLAC(object):
    def __init__(self):
        super(_Fma4BLAC,self).__init__()

    def Fma(self,s0Params,s1Params,s2Params,dParams,opts):

class FMA(ISA):
    def __init__(self,opts):
        super(FMA,self).__init__()
        self.name = "FMA"

        fp_m256d = { 'type' : '__m256d' }
        fp_m256d['arith'] = [ mm256FmaddPd , mm256FmaddsubPd , mm256FmsubPd , mm256FmsubaddPd , mm256FnmaddPd , mm256FnmsubPd ]
        fp_m256d['nublac'] = _Fma4BLAC()

        self.types = { 'fp' : { ('double', 4) : fp_m256d } } #only adding double precision data type
        #not sure what this line does
        #self.add_func_defs = [ asm256LoaduPd, asm256StoreuPd ]
