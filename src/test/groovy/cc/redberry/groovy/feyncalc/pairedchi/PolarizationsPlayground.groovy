/*
 * Redberry: symbolic tensor computations.
 *
 * Copyright (c) 2010-2015:
 *   Stanislav Poslavsky   <stvlpos@mail.ru>
 *   Bolotin Dmitriy       <bolotin.dmitriy@gmail.com>
 *
 * This file is part of Redberry.
 *
 * Redberry is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Redberry is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Redberry. If not, see <http://www.gnu.org/licenses/>.
 */
package cc.redberry.groovy.feyncalc.pairedchi

import cc.redberry.core.context.OutputFormat
import cc.redberry.groovy.Redberry
import org.junit.Test

import static cc.redberry.groovy.RedberryStatic.*

/**
 * Created by poslavsky on 04/06/15.
 */
class PolarizationsPlayground {
    private static def solverOptions = [ExternalSolver: [
            Solver: 'Mathematica',
            Path  : '/Applications/Mathematica.app/Contents/MacOS']
    ]

    static def gluonPolarizations(Setup stp) {
        use(Redberry) {
            def eps1 = 'eps1_a = c1 * k1_a + c2 * k2_a + c3 * p_a[bottom]'.t
            def eps2 = 'eps2_a = c4 * e_abcd * k1^b * k2^c * p^d[bottom]'.t
            def eqs = ['k1_a * eps1^a = 0'.t,
                       'k2_a * eps1^a = 0'.t,
                       'eps1_a * eps1^a = -1'.t,
                       'eps2_a * eps2^a = -1'.t,
                       'eps1_a * eps2^a = 0'.t]
            eqs = (eps1 & eps2 & ExpandAndEliminate & stp.leviSimplify &
                    ExpandAndEliminate & stp.mandelstam & stp.massesSubs) >> eqs

            def solutions = Reduce(eqs, ['c1', 'c2', 'c3', 'c4'].t, solverOptions)
            assert solutions.size() != 0

            def cfs = solutions[0]
            cfs = (stp.wolframFactorTr & stp.wolframFactorSqrtTr) >> cfs

            def epsPlus = 'eps_a[1] = (eps1_a + I * eps2_a)/2**(1/2)'.t
            def epsMinus = 'eps_a[-1] = (eps1_a - I * eps2_a)/2**(1/2)'.t
            return [epss: [eps1, eps2, epsPlus, epsMinus], coeffs: cfs]
        }
    }

    static def quarkoniaPolarizations(Setup stp, fl) {
        use(Redberry) {
            def var = fl[0] + 's'
            def eps1 = "eps1_a = ${var}1 * p_a[charm] + ${var}2 * p_a[bottom]".t
            def eps2 = "eps2_a = ${var}3 * p_a[charm] + ${var}4 * p_a[bottom] + ${var}5 * k1_a".t
            def eps0 = "eps0_a = ${var}6 * e_abcd * k1^b * p^c[charm] * p^d[bottom]".t

            if (!stp.projectCC) {
                def subs = 'p_a[charm] = p1_a[charm] + p2_a[charm]'.t.hold
                eps1 <<= subs; eps2 <<= subs; eps0 <<= subs;
            }

            def eq = ["p_a[$fl] * eps1^a = 0".t,
                      "p_a[$fl] * eps2^a = 0".t,
                      "p_a[$fl] * eps0^a = 0".t,
                      'eps1_a * eps2^a = 0'.t,
                      'eps1_a * eps1^a = -1'.t,
                      'eps2_a * eps2^a = -1'.t,
                      'eps0_a * eps0^a = -1'.t]
            eq = (eps1 & eps2 & eps0 & ExpandAndEliminate & stp.leviSimplify &
                    ExpandAndEliminate & stp.mandelstam & stp.massesSubs) >> eq
            def solutions = Reduce(eq, ["${var}1", "${var}2", "${var}3", "${var}4", "${var}5", "${var}6"].t, solverOptions)
            assert solutions.size() != 0

            def cfs = solutions[0]
            cfs = (stp.wolframFactorTr & stp.wolframFactorSqrtTr) >> cfs

            def epsPlus = "eps_a[$fl, 1] = (eps1_a + I * eps2_a)/2**(1/2)".t
            def epsMinus = "eps_a[$fl, -1] = (eps1_a - I * eps2_a)/2**(1/2)".t
            def epsZero = "eps_a[$fl, 0] = eps0_a".t

            //building tensor polarizations
            //def subs = Identity
            //for (def sub in ['eps1_a * eps1_b',
            //                 'eps2_a * eps2_b',
            //                 'eps1_a * eps2_b',
            //                 'eps1_a * eps0_b',
            //                 'eps2_a * eps0_b',
            //                 'eps0_a * eps0_b'].t) {
            //    def rhs = sub
            //    rhs <<= eps0 & eps1 & eps2 & cfs & stp.fullSimplifyE & stp.massesSubs & stp.wFactor
            //    subs &= sub.eq rhs
            //}
            def tEpss = []
            tEpss << "eps_ab[$fl, 2] = eps_a[$fl, 1] * eps_b[$fl, 1]".t
            tEpss << "eps_ab[$fl, 1] = (eps_a[$fl, 1]*eps_b[$fl, 0] + eps_a[$fl, 0]*eps_b[$fl, 1])/2**(1/2)".t
            tEpss << "eps_ab[$fl, 0] = 1/6**(1/2)*eps_a[$fl, 1]*eps_b[$fl, -1] + (2/3)**(1/2)*eps_a[$fl, 0]*eps_b[$fl, 0] + (1/6)**(1/2)*eps_a[$fl, -1]*eps_b[$fl, 1]".t
            tEpss << "eps_ab[$fl, -1] = (eps_a[$fl, -1] * eps_b[$fl, 0] + eps_a[$fl, 0]*eps_b[$fl, -1])/2**(1/2)".t
            tEpss << "eps_ab[$fl, -2] = eps_a[$fl, -1] * eps_b[$fl, -1]".t

            return [epss: [eps0, eps1, eps2, epsZero, epsPlus, epsMinus, *tEpss], coeffs: cfs]
        }
    }

    @Test
    public void testGluonPolarizations() throws Exception {
        Setup stp = new Setup(false, false, true)
        def gls = gluonPolarizations(stp)
        gls['coeffs'].each {
            println it.toString(OutputFormat.WolframMathematica)
        }
    }

    @Test
    public void testQuarkoniaPolarizations() throws Exception {
        Setup stp = new Setup(true, false, true)
        def coeffs = quarkoniaPolarizations(stp, 'charm')['coeffs']
        coeffs.each {
            println "def ${it[0]} = '$it'.t "
        }
    }

    static def coeffs = [:]

    private static def calcCoeffs(Setup stp, fl) {
        if (coeffs[fl] != null)
            return coeffs[fl]
        use(Redberry) {
            def var = fl[0] + 's'

            def eps1 = "eps1_a = ${var}1 * p_a[charm] + ${var}2 * p_a[bottom]".t
            def eps2 = "eps2_a = ${var}3 * p_a[charm] + ${var}4 * p_a[bottom] + ${var}5 * k1_a".t
            def eps0 = "eps0_a = ${var}6 * e_abcd * k1^b * p^c[charm] * p^d[bottom]".t

            if (!stp.projectCC) {
                def subs = 'p_a[charm] = p1_a[charm] + p2_a[charm]'.t.hold
                eps1 <<= subs; eps2 <<= subs; eps0 <<= subs;
            }

            def eq = ["p_a[$fl] * eps1^a = 0".t,
                      "p_a[$fl] * eps2^a = 0".t,
                      "p_a[$fl] * eps0^a = 0".t,
                      'eps1_a * eps2^a = 0'.t,
                      'eps1_a * eps1^a = -1'.t,
                      'eps2_a * eps2^a = -1'.t,
                      'eps0_a * eps0^a = -1'.t]
            eq = (eps1 & eps2 & eps0 & ExpandAndEliminate & stp.leviSimplify &
                    ExpandAndEliminate & stp.mandelstam & stp.massesSubs) >> eq
            def solutions = Reduce(eq, ["${var}1", "${var}2", "${var}3", "${var}4", "${var}5", "${var}6"].t, solverOptions)
            assert solutions.size() != 0

            def cfs = solutions[0]
            cfs = (stp.wolframFactorTr & stp.wolframFactorSqrtTr) >> cfs
            return (coeffs[fl] = cfs)
        }
    }

    @Test
    public void test1() throws Exception {
        Setup stp = new Setup(false, false, true)
        calcCoeffs(stp, 'bottom').each {
            println(it.toString(OutputFormat.WolframMathematica) + ';')
        }
    }

    private static def getVector(Setup stp, def fl, def xyz) {
        use(Redberry) {
            def var = fl[0] + 's', epss = [:]

            epss['x'] = "eps1_a = ${var}1 * p_a[charm] + ${var}2 * p_a[bottom]".t
            epss['y'] = "eps2_a = ${var}3 * p_a[charm] + ${var}4 * p_a[bottom] + ${var}5 * k1_a".t
            epss['z'] = "eps0_a = ${var}6 * e_abcd * k1^b * p^c[charm] * p^d[bottom]".t

            def coffs = calcCoeffs(stp, fl)
            def eps = epss[xyz][1] << coffs
            eps <<= Together
            return ['den': Denominator >> eps, 'num': Numerator >> eps]
        }
    }


    @Test
    public void testDefs() throws Exception {
        use(Redberry) {
            Setup stp = new Setup(false, false, true)
//            println getVector(stp, 'charm', 'x')
//            println getVector(stp, 'charm', 'y')
//            println getVector(stp, 'charm', 'z')

            println getVector(stp, 'bottom', 'x')
            println getVector(stp, 'bottom', 'y')
            println getVector(stp, 'bottom', 'z')
        }
    }

    @Test
    public void testDefsCC() throws Exception {
        use(Redberry) {
            println bind(['k_i[h]': 'p_j+f_j'])
        }
    }
//    @Test
//    public void testPolarizations1() throws Exception {
//        //gluon polarizations
//        use(Redberry) {
//            def stp = new Setup(true);
//            for (def g in [1, 2]) {
//                def sum = "eps${g}_a[1] * eps${g}_b[-1] + eps${g}_b[1] * eps${g}_a[-1]".t
//                sum <<= stp.polarisations & stp.fullSimplify & stp.massesSubs & 'u = 4*mc**2 + 4*mb**2 -s -t'.t
//                assert '2*s**(-1)*k2_{a}*k1_{b}+2*s**(-1)*k2_{b}*k1_{a}-g_{ab}'.t == (stp.wFactor >> sum)
//            }
//        }
//    }
//
//
//    @Test
//    public void testPolarizations1cc() throws Exception {
//        //gluon polarizations
//        use(Redberry) {
//            def stp = new Setup(false, true);
//            for (def g in [1, 2]) {
//                def sum = "eps${g}_a[1] * eps${g}_b[-1] + eps${g}_b[1] * eps${g}_a[-1]".t
//                sum <<= stp.polarisations & stp.fullSimplify & stp.massesSubs & stp.wFactor
//                assert '2*s**(-1)*k2_{a}*k1_{b}+2*s**(-1)*k2_{b}*k1_{a}-g_{ab}'.t == (stp.wFactor >> sum)
//            }
//        }
//    }
//
//    @Test
//    public void testPolarizations2() throws Exception {
//        //axial mesons polarizations
//        use(Redberry) {
//            def stp = new Setup(true);
//            for (def fl in ['charm', 'bottom']) {
//                def sum = "eps_a[$fl, 1] * eps_b[$fl, -1] + eps_a[$fl, 0] * eps_b[$fl, 0] + eps_b[$fl, 1] * eps_a[$fl, -1]".t
//                sum <<= stp.polarisations & stp.fullSimplify & stp.massesSubs & 'u = 4*mc**2 + 4*mb**2 -s -t'.t & stp.wFactor
//
//                def expected = "eps_a[h[$fl]] * eps_b[h[$fl]]".t
//                expected <<= stp.epsSum & stp.massesSubs
//                assert (sum - expected) == 0.t
//            }
//        }
//    }
//
//    @Test
//    public void testPolarizations2cc() throws Exception {
//        //axial mesons polarizations
//        use(Redberry) {
//            def stp = new Setup(false, true);
//            def fl = 'bottom'
//            def sum = "eps_a[$fl, 1] * eps_b[$fl, -1] + eps_a[$fl, 0] * eps_b[$fl, 0] + eps_b[$fl, 1] * eps_a[$fl, -1]".t
//            sum <<= stp.polarisations & stp.fullSimplify & stp.massesSubs & stp.wFactor
//
//            def expected = "eps_a[h[$fl]] * eps_b[h[$fl]]".t
//            expected <<= stp.epsSum & stp.massesSubs
//            assert (sum - expected) == 0.t
//        }
//    }
//
//
//    @Test
//    public void testPolarizations3() throws Exception {
//        //tensor mesons polarizations
//        use(Redberry) {
//            def stp = new Setup(true, true);
//            for (def fl in ['bottom']) {
//                def sum = new SumBuilder()
//                def pairs = [["eps_ab[$fl, 2]", "eps_cd[$fl, -2]"].t,
//                             ["eps_ab[$fl, 1]", "eps_cd[$fl, -1]"].t,
//                             ["eps_ab[$fl, 0]", "eps_cd[$fl, 0]"].t,
//                             ["eps_ab[$fl, -1]", "eps_cd[$fl, 1]"].t,
//                             ["eps_ab[$fl, -2]", "eps_cd[$fl, 2]"].t]
//
//                for (def pair in pairs) {
//                    def a = pair[0], b = pair[1]
//                    a <<= stp.polarisations; b <<= stp.polarisations;
//                    for (int i = 0; i < a.size(); ++i)
//                        for (int j = 0; j < b.size(); ++j) {
//                            Product p = a[i] * b[j]
//                            def ind = p.dataSubProduct
//                            ind <<= stp.fullSimplify & stp.massesSubs
//
//                            if (ind.class == Sum)
//                                sum << FastTensors.multiplySumElementsOnFactor(ind, p.indexlessSubProduct)
//                            else
//                                sum << ind * p.indexlessSubProduct
//                        }
//                }
//
//                //eps_ab * eps_cd
//                sum = sum.build()
//                println sum.size()
//
//                sum <<= 'u = 4*mc**2 + 4*mb**2 - s - t'.t & stp.wFactor & stp.wFactor
//                sum.each {
//                    println it
//                }
//            }
//        }
//    }
//
//    @Test
//    public void testPolarizations3сс() throws Exception {
//        //tensor mesons polarizations
//        use(Redberry) {
//            def stp = new Setup(false, true);
//            for (def fl in ['bottom']) {
//                def sum = new SumBuilder()
//                def pairs = [["eps_ab[$fl, 2]", "eps_cd[$fl, -2]"].t,
//                             ["eps_ab[$fl, 1]", "eps_cd[$fl, -1]"].t,
//                             ["eps_ab[$fl, 0]", "eps_cd[$fl, 0]"].t,
//                             ["eps_ab[$fl, -1]", "eps_cd[$fl, 1]"].t,
//                             ["eps_ab[$fl, -2]", "eps_cd[$fl, 2]"].t]
//
//                for (def pair in pairs) {
//                    def a = pair[0], b = pair[1]
//                    a <<= stp.polarisations; b <<= stp.polarisations;
//                    for (int i = 0; i < a.size(); ++i)
//                        for (int j = 0; j < b.size(); ++j) {
//                            Product p = a[i] * b[j]
//                            def ind = p.dataSubProduct
//                            ind <<= stp.fullSimplify & stp.massesSubs
//
//                            if (ind.class == Sum)
//                                sum << FastTensors.multiplySumElementsOnFactor(ind, p.indexlessSubProduct)
//                            else
//                                sum << ind * p.indexlessSubProduct
//                        }
//                }
//
//                //eps_ab * eps_cd
//                sum = sum.build()
//                println sum.size()
//
//                sum <<= stp.wFactor & stp.wFactor
//                sum.each {
//                    println it
//                }
//            }
//        }
//    }
//
//    @Test
//    public void testPolarizationsTensorSum() throws Exception {
//        //tensor mesons polarizations
//        use(Redberry) {
//            def stp = new Setup(true, true);
//            for (def fl in ['charm', 'bottom']) {
//                def pairs = [["eps_ab[$fl, 2]", "eps^ab[$fl, -2]"].t,
//                             ["eps_ab[$fl, 1]", "eps^ab[$fl, -1]"].t,
//                             ["eps_ab[$fl, 0]", "eps^ab[$fl, 0]"].t,
//                             ["eps_ab[$fl, -1]", "eps^ab[$fl, 1]"].t,
//                             ["eps_ab[$fl, -2]", "eps^ab[$fl, 2]"].t]
//
//                for (def pair in pairs) {
//                    def sum = new SumBuilder()
//                    def a = pair[0], b = pair[1]
//                    a <<= stp.polarisations; b <<= stp.polarisations;
//                    for (int i = 0; i < a.size(); ++i)
//                        for (int j = 0; j < b.size(); ++j) {
//                            Product p = a[i] * b[j]
//                            def ind = p.dataSubProduct
//                            ind <<= stp.fullSimplify & stp.massesSubs
//                            sum << p.indexlessSubProduct * ind
//                        }
//
//                    def res = sum.build()
//                    res <<= 'u = 4*mc**2 + 4*mb**2 -s -t'.t & stp.wolframFactorTr
//                    assert res == 1.t
//                }
//            }
//        }
//    }
//
//
//    @Test
//    public void testPolarizationsTensorSum2() throws Exception {
//        //tensor mesons polarizations
//        use(Redberry) {
//            def stp = new Setup(true, true);
//            for (def fl in ['charm', 'bottom']) {
//
//                for (def l1 in [-2, -1, 0, 1, 2])
//                    for (def l2 in [-2, -1, 0, 1, 2]) {
//                        def sum = new SumBuilder()
//                        def a = "eps_ab[$fl, $l1]".t, b = "eps^ab[$fl, $l2]".t
//                        a <<= stp.polarisations; b <<= stp.polarisations;
//                        for (int i = 0; i < a.size(); ++i)
//                            for (int j = 0; j < b.size(); ++j) {
//                                Product p = a[i] * b[j]
//                                def ind = p.dataSubProduct
//                                ind <<= stp.fullSimplify & stp.massesSubs
//                                sum << p.indexlessSubProduct * ind
//                            }
//                        def res = sum.build()
//                        res <<= 'u = 4*mc**2 + 4*mb**2 -s -t'.t & stp.wolframFactorTr
//                        if (l1 == -l2)
//                            assert res == 1.t
//                        else
//                            assert res == 0.t
//                    }
//            }
//        }
//    }
//
//    @Test
//    public void testPolarizationsTensorSum2cc() throws Exception {
//        //tensor mesons polarizations
//        use(Redberry) {
//            def stp = new Setup(false, true);
//            def fl = 'bottom'
//
//            for (def l1 in [-1, 0, 1, 2, -2])
//                for (def l2 in [-1, 0, 1, 2, -2]) {
//                    println "$l1 $l2"
//                    def sum = new SumBuilder()
//                    def a = "eps_ab[$fl, $l1]".t, b = "eps^ab[$fl, $l2]".t
//                    a <<= stp.polarisations; b <<= stp.polarisations;
//                    for (int i = 0; i < a.size(); ++i)
//                        for (int j = 0; j < b.size(); ++j) {
//                            Product p = a[i] * b[j]
//                            def ind = p.dataSubProduct
//                            ind <<= stp.fullSimplify & stp.massesSubs
//                            sum << p.indexlessSubProduct * ind
//                        }
//                    def res = sum.build()
//                    res <<= stp.wolframFactorTr
//                    if (l1 == -l2)
//                        assert res == 1.t
//                    else
//                        assert res == 0.t
//                }
//        }
//    }
//
//    @Test
//    public void testPolarizationsTensorSym() throws Exception {
//        //tensor mesons polarizations
//        use(Redberry) {
//            def stp = new Setup(true, true);
//            for (def fl in ['charm', 'bottom']) {
//
//                for (def l in [-2, -1, 0, 1, 2]) {
//                    def a = "eps_ab[$fl, $l] - eps_ba[$fl, $l]".t
//                    a <<= stp.polarisations & stp.wFactor
//                    assert a == 0.t
//
//                    a = "eps_ab[$fl, $l] * p^a[$fl]".t
//                    a <<= stp.polarisations & stp.fullSimplifyE & stp.massesSubs & stp.wFactor
//                    assert a == 0.t
//
//                    a = "eps_a^a[$fl, $l]".t
//                    a <<= stp.polarisations & stp.fullSimplifyE & stp.massesSubs
//                    a <<= 'u = 4*mc**2 + 4*mb**2 -s -t'.t & stp.wolframFactorTr
//                    assert a == 0.t
//                }
//            }
//        }
//    }
//
//    @Test
//    public void testPolarizationsTensorSymcc() throws Exception {
//        //tensor mesons polarizations
//        use(Redberry) {
//            def stp = new Setup(true, true);
//            def fl = 'bottom'
//            for (def l in [-2, -1, 0, 1, 2]) {
//                def a = "eps_ab[$fl, $l] - eps_ba[$fl, $l]".t
//                a <<= stp.polarisations & stp.wFactor
//                assert a == 0.t
//
//                a = "eps_ab[$fl, $l] * p^a[$fl]".t
//                a <<= stp.polarisations & stp.fullSimplifyE & stp.massesSubs & stp.wFactor
//                assert a == 0.t
//
//                a = "eps_a^a[$fl, $l]".t
//                a <<= stp.polarisations & stp.fullSimplifyE & stp.massesSubs
//                a <<= stp.wolframFactorTr
//                assert a == 0.t
//            }
//        }
//    }
}
