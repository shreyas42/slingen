'''
    Created on May 7, 2013
    
    @author: danieles
    '''

from islpy import Set, Map

from src.dsls.ll import Matrix, scalar_block
from src.physical import Array
from src.binding import getReference, ExplicitPhysicalReference

from src.irbase import icode
# from src.physical import *
# from src.binding import *
# import src.irbase as irbase
# from irbase import *

class Allocator(object):
    def __init__(self):
        super(Allocator, self).__init__()
    
    def declare(self, m, subM, opts):
        if subM is not None and icode.bindingTable.isBound(subM):
            subPhys = icode.bindingTable.getPhysicalLayout(subM)
            icode.bindingTable.addBinding(m, subPhys)
        else:
            #changing this for now : MARKER
            is_complex = False
            if m.get_field() == 'complex':
                is_complex = True
            if is_complex:
                outPhys = Array(m.name, m.size, opts, safelyScalarize=opts['scarep'], field = 'BlkInterLeaved')
            else:
                outPhys = Array(m.name, m.size, opts, safelyScalarize=opts['scarep'])
            if icode.bindingTable.addBinding(m, outPhys):
                icode.declare += [outPhys]
    
    def reshape(self, m, M, N, bounds, opts):
        # With Heterogeneous layouts, the size of certain base cases may depend upon some of the indices.
        # In MMM, the "accumulation" matrix's PhysLayout is one such a case.
        currPhys = icode.getPhysicalLayout(m)
        if isinstance(currPhys, Array):
            if any(map(lambda idx: currPhys.hasIdx(idx), bounds)):
                currPhys.subs(bounds)
            elif currPhys.size < M*N:
                newPhys = Array(m.name, (M,N), opts, safelyScalarize=opts['scarep'])
                icode.replacePhysicalLayout(currPhys, newPhys)
    
    def extractParams(self, m, bounds, opts, bReshape=True, scaRowExt=True, scaColExt=False, rowExt=False, colExt=False):
        '''
            Return index mapping functions, dimensions, nu-niceness and pointer to nu-compliant matrix
            '''
        raise NotImplementedError('Should have implemented this!')


class NuAllocator(Allocator):
    def __init__(self):
        super(NuAllocator, self).__init__()
    
    def extractParams(self, m, bounds, opts, bReshape=False, scaRowExt=True, scaColExt=False, rowExt=False, colExt=False, bcast=False, subNuM=None):
        '''
            Return index mapping functions, dimensions, nu-niceness and pointer to nu-compliant matrix
            '''
        mL,mR = m.fL.subs(bounds), m.fR.subs(bounds)
        M, N = m.size[0].subs(bounds), m.size[1].subs(bounds)
        
        if bReshape:
            self.reshape(m, M, N, bounds, opts)

        ref = getReference(icode, m)
        #change added here
        physLayout = icode.bindingTable.getPhysicalLayout(m)
        is_exp_phys_ref = isinstance(ref, ExplicitPhysicalReference)

        is_complex_layout = False
        if physLayout.getField() in ['BlkInterLeaved','Split']:
            is_complex_layout = True

        isCompact, isCorner = None, None
        if is_exp_phys_ref:
            if is_complex_layout:
            #changed the two lines below , this needs to actually reflect what kind of Reference it is
                isCompact = (ref.getLinIdx([mL.of(M-1), mR.of(N-1) , 1]) - ref.getLinIdx([mL.of(0), mR.of(0) , 0])) == physLayout.getSize()
                isCorner = ref.isCorner([mL.of(M-1), mR.of(N-1) , 1])
            else:
                isCompact = (ref.getLinIdx([mL.of(M-1), mR.of(N-1)]) - ref.getLinIdx([mL.of(0), mR.of(0)])) == physLayout.getSize()
                isCorner = ref.isCorner([mL.of(M-1), mR.of(N-1)])


        access = m.genAccess()
        
        nu = opts['nu']
#         if not opts['vectorize'] or M > nu or N > nu:
#             return { 'mL': mL, 'mR': mR, 'M': M, 'N': N, 'nuable': False, 'struct': m.genStruct, 'access': access, 'nuM': None, 'm': m, 'nu': 0, 'compact': isCompact, 'corner': isCorner, 'bounds': bounds }
        if opts['useintrinsics'] and (M > nu or N > nu):
            return { 'mL': mL, 'mR': mR, 'M': M, 'N': N, 'nuable': False, 'struct': m.genStruct, 'access': access, 'nuM': None, 'm': m, 'nu': nu, 'compact': isCompact, 'corner': isCorner, 'bounds': bounds }

#         isnuGeMat = Matrix in m.genStruct and m.genStruct[Matrix] == Set("{[i,j]: 0<=i,j<%d}"%nu)
#         isnuGeMat = isnuGeMat and access.intersect(Map("{[i,j]->[i,j]}")) == access
        isnuGeMat = Matrix.testGeneral(m.genStruct, access, nu, nu)

        if (opts['useintrinsics'] or not is_exp_phys_ref) and M == 1 and N == 1:
            nuMM, nuMN = nu if scaColExt else 1, nu if scaRowExt else 1
            #change here : IMPORTANT CHANGE
            if is_complex_layout:
                nuM = Matrix('', scalar_block('complex'), (nuMM,nuMN))
            else:
                nuM = Matrix('', scalar_block(), (nuMM,nuMN))

            nuML, nuMR = nuM.fL, nuM.fR
            self.declare(nuM, subNuM, opts)
        elif is_exp_phys_ref and (not opts['useintrinsics'] or (M == nu and N == nu and isnuGeMat) or (M == 1 and N == nu) or (M == nu and N == 1 and isCompact)):
            nuMM, nuMN = M, N
            nuML, nuMR = mL, mR
            nuM = m
        elif M == 1 and N < nu:
            nuMM, nuMN = nu if colExt else 1, nu
            if is_complex_layout:
                nuM = Matrix('', scalar_block('complex'), (nuMM,nuMN))
            else:
                nuM = Matrix('', scalar_block(), (nuMM,nuMN))

            nuML, nuMR = nuM.fL, nuM.fR
            self.declare(nuM, subNuM, opts)
        elif M <= nu and N == 1:
            nuMM, nuMN = nu, nu if rowExt else 1
            if is_complex_layout:
                nuM = Matrix('', scalar_block('complex'), (nuMM,nuMN))
            else:
                nuM = Matrix('', scalar_block(), (nuMM,nuMN))

            nuML, nuMR = nuM.fL, nuM.fR
            self.declare(nuM, subNuM, opts)
        else:
            nuMM, nuMN = nu, nu
            if is_complex_layout:
                nuM = Matrix('', scalar_block('complex'), (nu,nu))
            else:
                nuM = Matrix('', scalar_block(), (nu,nu))

            nuML, nuMR = nuM.fL, nuM.fR
            self.declare(nuM, subNuM, opts)
        
        return { 'mL': mL, 'mR': mR, 'M': M, 'N': N, 'nuML': nuML, 'nuMR': nuMR, 'nuMM': nuMM, 'nuMN': nuMN, 'nuable': True, 'nuM': nuM, 'm': m, 'struct': m.genStruct, 'access': access, 'nu': nu, 'compact': isCompact, 'corner': isCorner, 'bounds': bounds, 'bcast': bcast }


class NeonAllocator(Allocator):
    def __init__(self):
        super(NeonAllocator, self).__init__()
    
    def extractParams(self, m, bounds, opts, bReshape=False, scaRowExt=True, scaColExt=False, rowExt=False, colExt=False):
        '''
            Return index mapping functions, dimensions, nu-niceness and pointer to nu-compliant matrix.
            
            Attention: Currently scaRowExt, scaColExt, rowExt, colExt are not taken into account.
            '''
        mL,mR = m.fL.subs(bounds), m.fR.subs(bounds)
        M, N = m.size[0].subs(bounds), m.size[1].subs(bounds)
        
        if bReshape:
            self.reshape(m, M, N, bounds, opts)
        
        if not opts['useintrinsics'] or M > 4 or N > 4:
            return { 'mL': mL, 'mR': mR, 'M': M, 'N': N, 'nuable': False, 'nuM': None, 'm': m, 'nu': 0, 'bounds': bounds }
        
        ref = getReference(icode, m)
        
        isCompact = (ref.getLinIdx([mL.of(M-1), mR.of(N-1)]) - ref.getLinIdx([mL.of(0), mR.of(0)])) == (M*N - 1)
        isCorner = ref.isCorner([mL.of(M-1), mR.of(N-1)])
        if N == 2 and (M == 1 or not isCompact) or M == 2 and N == 1 and isCompact or N == 4 or M == 4 and N == 1 and isCompact:
            nuMM, nuMN = M, N
            nuML, nuMR = mL, mR
            nuM = m
        elif M == 1:
            if N == 1:
                nuMM, nuMN = 2 if scaColExt else 1, 2 if scaRowExt else 1
                nuM = Matrix('', scalar_block(), (nuMM,nuMN))
                nuML, nuMR = nuM.fL, nuM.fR
                self.declare(nuM, opts)
            elif N == 3:
                nuMM, nuMN = 1, 4
                nuM = Matrix('', scalar_block(), (nuMM,nuMN))
                nuML, nuMR = nuM.fL, nuM.fR
                self.declare(nuM, opts)
        elif M == 2:
            if N == 1 and not isCompact:
                nuMM, nuMN = 1, 2
                nuM = Matrix('', scalar_block(), (nuMM,nuMN))
                nuML, nuMR = nuM.fL, nuM.fR
                self.declare(nuM, opts)
            elif N == 2 and isCompact:
                nuMM, nuMN = 2, 2
                nuM = Matrix('', scalar_block(), (nuMM,nuMN))
                nuML, nuMR = nuM.fL, nuM.fR
                self.declare(nuM, opts)
            elif N == 3:
                nuMM, nuMN = 2, 4
                nuM = Matrix('', scalar_block(), (nuMM,nuMN))
                nuML, nuMR = nuM.fL, nuM.fR
                self.declare(nuM, opts)
        elif M == 3:
            if N == 1:
                nuMM, nuMN = 1, 4
                nuM = Matrix('', scalar_block(), (nuMM,nuMN))
                nuML, nuMR = nuM.fL, nuM.fR
                self.declare(nuM, opts)
            elif N == 2 and isCompact:
                nuMM, nuMN = 3, 2
                nuM = Matrix('', scalar_block(), (nuMM,nuMN))
                nuML, nuMR = nuM.fL, nuM.fR
                self.declare(nuM, opts)
            elif N == 3:
                nuMM, nuMN = 3, 4
                nuM = Matrix('', scalar_block(), (nuMM,nuMN))
                nuML, nuMR = nuM.fL, nuM.fR
                self.declare(nuM, opts)
        elif M == 4:
            if N == 1 and not isCompact:
                nuMM, nuMN = 1, 4
                nuM = Matrix('', scalar_block(), (nuMM,nuMN))
                nuML, nuMR = nuM.fL, nuM.fR
                self.declare(nuM, opts)
            elif N == 2 and isCompact:
                nuMM, nuMN = 4, 2
                nuM = Matrix('', scalar_block(), (nuMM,nuMN))
                nuML, nuMR = nuM.fL, nuM.fR
                self.declare(nuM, opts)
            elif N == 3:
                nuMM, nuMN = 4, 4
                nuM = Matrix('', scalar_block(), (nuMM,nuMN))
                nuML, nuMR = nuM.fL, nuM.fR
                self.declare(nuM, opts)
        else:
            raise Exception('Error in allocator: M=%d, N=%d, isCompact=%s' % (M, N, str(isCompact)))
        return { 'mL': mL, 'mR': mR, 'M': M, 'N': N, 'nuML': nuML, 'nuMR': nuMR, 'nuMM': nuMM, 'nuMN': nuMN, 'nuable': True, 'nuM': nuM, 'm': m, 'nu': 4, 'compact': isCompact, 'corner': isCorner, 'bounds': bounds }
