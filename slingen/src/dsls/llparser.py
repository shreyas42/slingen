# 2017-08-18 10:44:14.284057
#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# CAVEAT UTILITOR
# This file was automatically generated by Grako.
#    https://bitbucket.org/apalala/grako/
# Any changes you make to it will be overwritten the
# next time the file is generated.
#

from __future__ import print_function, division, absolute_import, unicode_literals
from grako.parsing import *  # noqa
from grako.exceptions import *  # noqa


__version__ = '17.230.08.44.14'


class llParser(Parser):
    def __init__(self, whitespace=None, **kwargs):
        super(llParser, self).__init__(whitespace=whitespace, **kwargs)

    @rule_def
    def _program_(self):

        def block0():
            self._declaration_()
            self.ast.add_list('decl', self.last_node)
        self._closure(block0)

        def block2():
            self._statement_()
            self.ast.add_list('stmt', self.last_node)
        self._closure(block2)
        self._check_eof()

    @rule_def
    def _declaration_(self):
        with self._choice():
            with self._option():
                with self._group():
                    self._token('Scalar')
                    self.ast['vartype'] = self.last_node
                    self._id_()
                    self.ast['name'] = self.last_node
                    self._token('<')
                    self._io_()
                    self.ast['iotype'] = self.last_node
                    self._token('>')
                    self._token(';')
            with self._option():
                with self._group():
                    self._token('Vector')
                    self.ast['vartype'] = self.last_node
                    self._id_()
                    self.ast['name'] = self.last_node
                    self._dim_vector_()
                    self.ast['dims'] = self.last_node
                    self._token('<')
                    self._io_()
                    self.ast['iotype'] = self.last_node
                    with self._optional():
                        self._token(',')
                        self._ow_()
                        self.ast['ow'] = self.last_node
                    self._token('>')
                    self._token(';')
            with self._option():
                with self._group():
                    self._token('Matrix')
                    self.ast['vartype'] = self.last_node
                    self._id_()
                    self.ast['name'] = self.last_node
                    self._dim_matrix_()
                    self.ast['dims'] = self.last_node
                    self._token('<')
                    self._io_()
                    self.ast['iotype'] = self.last_node

                    def block12():
                        self._token(',')
                        self._prop_()
                        self.ast.add_list('props', self.last_node)
                    self._closure(block12)
                    with self._optional():
                        self._token(',')
                        self._ow_()
                        self.ast['ow'] = self.last_node
                    self._token('>')
                    self._token(';')
            self._error('no available options')

    @rule_def
    def _dim_vector_(self):
        self._token('(')
        with self._optional():
            self._id_()
            self.ast.add_list('id', self.last_node)
            self._token('=')
        self._numexpr_()
        self.ast.add_list('val', self.last_node)
        self._token(')')

    @rule_def
    def _dim_matrix_(self):
        self._token('(')
        with self._optional():
            self._id_()
            self.ast.add_list('id', self.last_node)
            self._token('=')
        self._numexpr_()
        self.ast.add_list('val', self.last_node)
        self._token(',')
        with self._optional():
            self._id_()
            self.ast.add_list('id', self.last_node)
            self._token('=')
        self._numexpr_()
        self.ast.add_list('val', self.last_node)
        self._token(')')

    @rule_def
    def _io_(self):
        with self._choice():
            with self._option():
                self._token('Input')
            with self._option():
                self._token('Output')
            with self._option():
                self._token('InOut')
            with self._option():
                self._token('tInput')
            with self._option():
                self._token('tOutput')
            with self._option():
                self._token('tInOut')
            self._error('expecting one of: Input tInput tInOut tOutput InOut Output')

    @rule_def
    def _prop_(self):
        with self._choice():
            with self._option():
                self._token('Square')
                self.ast['@'] = self.last_node
            with self._option():
                self._token('Rectangular')
                self.ast['@'] = self.last_node
            with self._option():
                self._token('Diagonal')
                self.ast['@'] = self.last_node
            with self._option():
                self._token('LowerTriangular')
                self.ast['@'] = self.last_node
            with self._option():
                self._token('UpperTriangular')
                self.ast['@'] = self.last_node
            with self._option():
                self._token('UnitDiagonal')
                self.ast['@'] = self.last_node
            with self._option():
                self._token('ImplicitUnitDiagonal')
                self.ast['@'] = self.last_node
            with self._option():
                self._token('Symmetric')
                self.ast['@'] = self.last_node
            with self._option():
                self._token('Non-singular')
                self.ast['@'] = self.last_node
            with self._option():
                self._token('SPD')
                self.ast['@'] = self.last_node
            with self._option():
                self._token('LowerStorage')
                self.ast['@'] = self.last_node
            with self._option():
                self._token('UpperStorage')
                self.ast['@'] = self.last_node
            self._error('expecting one of: Rectangular UpperStorage SPD Square Non-singular Symmetric LowerTriangular UpperTriangular UnitDiagonal ImplicitUnitDiagonal Diagonal LowerStorage')

    @rule_def
    def _ow_(self):
        self._token('overwrites')
        self._token('(')
        self._id_()
        self.ast['@'] = self.last_node
        self._token(')')

    @rule_def
    def _statement_(self):
        with self._group():
            with self._choice():
                with self._option():
                    self._llfor_()
                with self._option():
                    self._llif_()
                with self._option():
                    self._equation_()
                self._error('no available options')
        self._token(';')

    @rule_def
    def _llfor_(self):
        self._token('For')
        self._cut()
        self._preprocs_()
        self._looptop_()
        self.ast['looptop'] = self.last_node
        self._token('{')

        def block1():
            self._statement_()
            self.ast.add_list('body', self.last_node)
        self._positive_closure(block1)

        self._token('}')

    @rule_def
    def _looptop_(self):
        self._token('[')
        self._id_()
        self.ast['idx'] = self.last_node
        self._token(';')
        self._numexpr_()
        self.ast['lb'] = self.last_node
        self._token(';')
        self._numexpr_()
        self.ast['ub'] = self.last_node
        self._token(';')
        self._numexpr_()
        self.ast['s'] = self.last_node
        self._token(']')

    @rule_def
    def _llif_(self):
        self._token('If')
        self._cut()
        self._preprocs_()
        self._guard_()
        self.ast['guard'] = self.last_node
        self._token('{')

        def block1():
            self._statement_()
            self.ast.add_list('then', self.last_node)
        self._positive_closure(block1)

        self._token('}')

    @rule_def
    def _guard_(self):
        self._token('[')
        self._condexpr_()
        self.ast['ce'] = self.last_node
        self._token(']')

    @rule_def
    def _equation_(self):
        with self._optional():
            self._eqann_()
            self.ast['eqann'] = self.last_node
        self._lhs_()
        self.ast['lhs'] = self.last_node
        self._token('=')
        self._rhs_()
        self.ast['rhs'] = self.last_node

    @rule_def
    def _eqann_(self):
        self._token('#')
        self._eqanntxt_()
        self.ast.add_list('@', self.last_node)

        def block1():
            self._token(';')
            self._eqanntxt_()
            self.ast.add_list('@', self.last_node)
        self._closure(block1)
        self._token('#')

    @rule_def
    def _eqanntxt_(self):
        self._pattern(r'[^;#]*')

    @rule_def
    def _lhs_(self):
        with self._choice():
            with self._option():
                self._lexpr_()
                self.ast['@'] = self.last_node
            with self._option():
                self._token('[')
                self._lexpr_()
                self.ast.add_list('multiout', self.last_node)

                def block2():
                    self._token(',')
                    self._lexpr_()
                    self.ast.add_list('multiout', self.last_node)
                self._closure(block2)
                self._token(']')
            self._error('no available options')

    @rule_def
    def _lexpr_(self):
        self._lterm_()
        self.ast.add_list('term', self.last_node)

        def block1():
            self._add_()
            self.ast.add_list('op', self.last_node)
            self._lterm_()
            self.ast.add_list('term', self.last_node)
        self._closure(block1)

    @rule_def
    def _lterm_(self):
        self._planefactor_()
        self.ast.add_list('factor', self.last_node)

        def block1():
            self._md_()
            self.ast.add_list('fop', self.last_node)
            self._cut()
            self._planefactor_()
            self.ast.add_list('factor', self.last_node)
        self._closure(block1)

    @rule_def
    def _rhs_(self):
        self._expr_()
        self.ast['@'] = self.last_node

    @rule_def
    def _expr_(self):
        self._term_()
        self.ast.add_list('term', self.last_node)

        def block1():
            self._add_()
            self.ast.add_list('op', self.last_node)
            self._term_()
            self.ast.add_list('term', self.last_node)
        self._closure(block1)

    @rule_def
    def _term_(self):
        self._factor_()
        self.ast.add_list('factor', self.last_node)

        def block1():
            self._md_()
            self.ast.add_list('fop', self.last_node)
            self._cut()
            self._factor_()
            self.ast.add_list('factor', self.last_node)
        self._closure(block1)

    @rule_def
    def _factor_(self):
        with self._choice():
            with self._option():
                self._scatter_()
                self.ast['@'] = self.last_node
            with self._option():
                self._scatteracc_()
                self.ast['@'] = self.last_node
            with self._option():
                self._planefactor_()
                self.ast['@'] = self.last_node
            self._error('no available options')

    @rule_def
    def _scatter_(self):
        self._token('[')
        self._imf_()
        self.ast['simf'] = self.last_node
        self._token(',')
        self._imf_()
        self.ast['simf'] = self.last_node
        self._token(']')
        self._factor_()
        self.ast['factor'] = self.last_node

    @rule_def
    def _scatteracc_(self):
        self._token('$[')
        self._imf_()
        self.ast['saimf'] = self.last_node
        self._token(',')
        self._imf_()
        self.ast['saimf'] = self.last_node
        self._token(']')
        self._factor_()
        self.ast['factor'] = self.last_node

    @rule_def
    def _planefactor_(self):
        with self._optional():
            self._add_()
            self.ast['sign'] = self.last_node
        with self._group():
            with self._choice():
                with self._option():
                    self._token('inv')
                    self._token('(')
                    self._expr_()
                    self.ast['inv'] = self.last_node
                    self._token(')')
                with self._option():
                    self._token('trans')
                    self._token('(')
                    self._expr_()
                    self.ast['trans'] = self.last_node
                    self._token(')')
                with self._option():
                    self._token('sqrt')
                    self._token('(')
                    self._expr_()
                    self.ast['sqrt'] = self.last_node
                    self._token(')')
                with self._option():
                    self._funcall_()
                    self.ast['func'] = self.last_node
                with self._option():
                    self._token('(')
                    self._expr_()
                    self.ast['par'] = self.last_node
                    self._token(')')
                with self._option():
                    self._id_()
                    self.ast['id'] = self.last_node
                with self._option():
                    self._constnum_()
                    self.ast['const'] = self.last_node
                self._error('no available options')
        self._preprocg_()
        self._gather_()
        self.ast['gather'] = self.last_node

    @rule_def
    def _funcall_(self):
        self._id_()
        self.ast['name'] = self.last_node
        self._token('(')
        self._numexpr_()
        self.ast['m'] = self.last_node
        self._token(',')
        self._numexpr_()
        self.ast['n'] = self.last_node
        self._token(';')
        self._expr_()
        self.ast.add_list('params', self.last_node)

        def block4():
            self._token(',')
            self._cut()
            self._expr_()
            self.ast.add_list('params', self.last_node)
        self._closure(block4)
        self._token(')')

    @rule_def
    def _gather_(self):
        with self._choice():
            with self._option():
                self._token('[')
                self._imf_()
                self._token(',')
                self._imf_()
                self._token(']')
                with self._optional():
                    self._ann_()
                    self.ast['ann'] = self.last_node
                self._gather_()
            with self._option():
                pass
            self._error('no available options')

    @rule_def
    def _ann_(self):
        self._token('#')
        self._txt_()
        self.ast.add_list('@', self.last_node)

        def block1():
            self._token(',')
            self._txt_()
            self.ast.add_list('@', self.last_node)
        self._closure(block1)
        self._token('#')

    @rule_def
    def _preprocg_(self):
        pass

    @rule_def
    def _preprocs_(self):
        pass

    @rule_def
    def _imf_(self):
        with self._choice():
            with self._option():
                self._iimf_()
                self.ast['@'] = self.last_node
            with self._option():
                self._genimf_()
                self.ast['@'] = self.last_node
            with self._option():
                self._himf_()
                self.ast['@'] = self.last_node
            self._error('no available options')

    @rule_def
    def _genimf_(self):
        self._token('f(')
        self._cut()
        self._numexpr_()
        self.ast['params'] = self.last_node
        self._token(',')
        self._numexpr_()
        self.ast['params'] = self.last_node
        self._token(',')
        self._numexpr_()
        self.ast['params'] = self.last_node
        self._token(',')
        self._numexpr_()
        self.ast['params'] = self.last_node
        self._token(')')

    @rule_def
    def _himf_(self):
        self._token('h(')
        self._cut()
        self._numexpr_()
        self.ast['params'] = self.last_node
        self._token(',')
        self._numexpr_()
        self.ast['params'] = self.last_node
        self._token(',')
        self._numexpr_()
        self.ast['params'] = self.last_node
        with self._optional():
            self._token(',')
            self._numexpr_()
            self.ast['params'] = self.last_node
        self._token(')')

    @rule_def
    def _iimf_(self):
        self._token('fI(')
        self._cut()
        self._numexpr_()
        self.ast['params'] = self.last_node
        self._token(')')

    @rule_def
    def _condexpr_(self):
        self._condterm_()
        self.ast.add_list('condterm', self.last_node)

        def block1():
            self._token('||')
            self._cut()
            self._condterm_()
            self.ast.add_list('condterm', self.last_node)
        self._closure(block1)

    @rule_def
    def _condterm_(self):
        self._condfactor_()
        self.ast.add_list('condfactor', self.last_node)

        def block1():
            self._token('&&')
            self._cut()
            self._condfactor_()
            self.ast.add_list('condfactor', self.last_node)
        self._closure(block1)

    @rule_def
    def _condfactor_(self):
        with self._optional():
            self._token('(')
        self._numexpr_()
        self.ast['condl'] = self.last_node
        self._condsym_()
        self.ast['condsym'] = self.last_node
        self._numexpr_()
        self.ast['condr'] = self.last_node
        with self._optional():
            self._token(')')

    @rule_def
    def _numexpr_(self):
        with self._optional():
            self._add_()
            self.ast['sign'] = self.last_node
        self._numterm_()
        self.ast.add_list('numterm', self.last_node)

        def block2():
            self._add_()
            self.ast.add_list('op', self.last_node)
            self._cut()
            self._numterm_()
            self.ast.add_list('numterm', self.last_node)
        self._closure(block2)

    @rule_def
    def _numterm_(self):
        self._numfactor_()
        self.ast.add_list('numfactor', self.last_node)

        def block1():
            self._mul_()
            self.ast.add_list('op', self.last_node)
            self._cut()
            self._numfactor_()
            self.ast.add_list('numfactor', self.last_node)
        self._closure(block1)

    @rule_def
    def _numfactor_(self):
        with self._choice():
            with self._option():
                self._mod_()
            with self._option():
                self._floord_()
            with self._option():
                self._ceild_()
            with self._option():
                self._min_()
            with self._option():
                self._max_()
            with self._option():
                self._id_()
                self.ast['id'] = self.last_node
            with self._option():
                self._zint_()
                self.ast['const'] = self.last_node
            with self._option():
                self._numsubexpr_()
            self._error('no available options')

    @rule_def
    def _mod_(self):
        with self._group():
            with self._choice():
                with self._option():
                    self._token('mod(')
                with self._option():
                    self._token('Mod(')
                self._error('expecting one of: mod( Mod(')
        self._cut()
        self._numexpr_()
        self.ast['modl'] = self.last_node
        self._token(',')
        self._numexpr_()
        self.ast['modr'] = self.last_node
        self._token(')')

    @rule_def
    def _floord_(self):
        self._token('floord(')
        self._cut()
        self._numexpr_()
        self.ast['fnum'] = self.last_node
        self._token(',')
        self._int_()
        self.ast['fden'] = self.last_node
        self._token(')')

    @rule_def
    def _numsubexpr_(self):
        self._token('(')
        self._cut()
        self._numexpr_()
        self.ast['numpar'] = self.last_node
        self._token(')')

    @rule_def
    def _ceild_(self):
        self._token('ceild(')
        self._cut()
        self._numexpr_()
        self.ast['cnum'] = self.last_node
        self._token(',')
        self._int_()
        self.ast['cden'] = self.last_node
        self._token(')')

    @rule_def
    def _min_(self):
        with self._group():
            with self._choice():
                with self._option():
                    self._token('min(')
                with self._option():
                    self._token('Min(')
                self._error('expecting one of: min( Min(')
        self._cut()
        self._numexpr_()
        self.ast['minl'] = self.last_node
        self._token(',')
        self._numexpr_()
        self.ast['minr'] = self.last_node
        self._token(')')

    @rule_def
    def _max_(self):
        with self._group():
            with self._choice():
                with self._option():
                    self._token('max(')
                with self._option():
                    self._token('Max(')
                self._error('expecting one of: max( Max(')
        self._cut()
        self._numexpr_()
        self.ast['maxl'] = self.last_node
        self._token(',')
        self._numexpr_()
        self.ast['maxr'] = self.last_node
        self._token(')')

    @rule_def
    def _condsym_(self):
        with self._choice():
            with self._option():
                self._token('==')
            with self._option():
                self._token('>=')
            with self._option():
                self._token('>')
            with self._option():
                self._token('<=')
            with self._option():
                self._token('<')
            with self._option():
                self._token('!=')
            self._error('expecting one of: > >= <= != < ==')

    @rule_def
    def _add_(self):
        self._pattern(r'[\+-]')

    @rule_def
    def _mul_(self):
        self._pattern(r'[\*/]')

    @rule_def
    def _md_(self):
        self._pattern(r'[\\\*/]')

    @rule_def
    def _txt_(self):
        self._pattern(r'[^,#]*')

    @rule_def
    def _constnum_(self):
        self._pattern(r'[0-9]+')

    @rule_def
    def _zint_(self):
        with self._choice():
            with self._option():
                self._token('0')
            with self._option():
                self._int_()
            self._error('expecting one of: 0')

    @rule_def
    def _int_(self):
        with self._choice():
            with self._option():
                self._pattern(r'[1-9][0-9]*')
            with self._option():
                self._pattern(r'[A-Za-z]')
            self._error('expecting one of: [A-Za-z] [1-9][0-9]*')

    @rule_def
    def _id_(self):
        self._pattern(r'[A-Za-z][A-Za-z0-9_]*')


class llSemanticParser(CheckSemanticsMixin, llParser):
    pass


class llSemantics(object):
    def program(self, ast):
        return ast

    def declaration(self, ast):
        return ast

    def dim_vector(self, ast):
        return ast

    def dim_matrix(self, ast):
        return ast

    def io(self, ast):
        return ast

    def prop(self, ast):
        return ast

    def ow(self, ast):
        return ast

    def statement(self, ast):
        return ast

    def llfor(self, ast):
        return ast

    def looptop(self, ast):
        return ast

    def llif(self, ast):
        return ast

    def guard(self, ast):
        return ast

    def equation(self, ast):
        return ast

    def eqann(self, ast):
        return ast

    def eqanntxt(self, ast):
        return ast

    def lhs(self, ast):
        return ast

    def lexpr(self, ast):
        return ast

    def lterm(self, ast):
        return ast

    def rhs(self, ast):
        return ast

    def expr(self, ast):
        return ast

    def term(self, ast):
        return ast

    def factor(self, ast):
        return ast

    def scatter(self, ast):
        return ast

    def scatteracc(self, ast):
        return ast

    def planefactor(self, ast):
        return ast

    def funcall(self, ast):
        return ast

    def gather(self, ast):
        return ast

    def ann(self, ast):
        return ast

    def preprocg(self, ast):
        return ast

    def preprocs(self, ast):
        return ast

    def imf(self, ast):
        return ast

    def genimf(self, ast):
        return ast

    def himf(self, ast):
        return ast

    def iimf(self, ast):
        return ast

    def condexpr(self, ast):
        return ast

    def condterm(self, ast):
        return ast

    def condfactor(self, ast):
        return ast

    def numexpr(self, ast):
        return ast

    def numterm(self, ast):
        return ast

    def numfactor(self, ast):
        return ast

    def mod(self, ast):
        return ast

    def floord(self, ast):
        return ast

    def numsubexpr(self, ast):
        return ast

    def ceild(self, ast):
        return ast

    def min(self, ast):
        return ast

    def max(self, ast):
        return ast

    def condsym(self, ast):
        return ast

    def add(self, ast):
        return ast

    def mul(self, ast):
        return ast

    def md(self, ast):
        return ast

    def txt(self, ast):
        return ast

    def constnum(self, ast):
        return ast

    def zint(self, ast):
        return ast

    def int(self, ast):
        return ast

    def id(self, ast):
        return ast


def main(filename, startrule, trace=False, whitespace=None):
    import json
    with open(filename) as f:
        text = f.read()
    parser = llParser(parseinfo=False)
    ast = parser.parse(
        text,
        startrule,
        filename=filename,
        trace=trace,
        whitespace=whitespace)
    print('AST:')
    print(ast)
    print()
    print('JSON:')
    print(json.dumps(ast, indent=2))
    print()

if __name__ == '__main__':
    import argparse
    import string
    import sys

    class ListRules(argparse.Action):
        def __call__(self, parser, namespace, values, option_string):
            print('Rules:')
            for r in llParser.rule_list():
                print(r)
            print()
            sys.exit(0)

    parser = argparse.ArgumentParser(description="Simple parser for ll.")
    parser.add_argument('-l', '--list', action=ListRules, nargs=0,
                        help="list all rules and exit")
    parser.add_argument('-t', '--trace', action='store_true',
                        help="output trace information")
    parser.add_argument('-w', '--whitespace', type=str, default=string.whitespace,
                        help="whitespace specification")
    parser.add_argument('file', metavar="FILE", help="the input file to parse")
    parser.add_argument('startrule', metavar="STARTRULE",
                        help="the start rule for parsing")
    args = parser.parse_args()

    main(args.file, args.startrule, trace=args.trace, whitespace=args.whitespace)