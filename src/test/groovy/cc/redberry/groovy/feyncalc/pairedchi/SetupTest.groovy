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

import cc.redberry.core.context.CC
import cc.redberry.core.tensor.SumBuilder
import cc.redberry.core.transformations.Transformation
import cc.redberry.core.utils.TensorUtils
import cc.redberry.groovy.Redberry
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics
import org.junit.Before
import org.junit.Ignore
import org.junit.Test

import static cc.redberry.core.indices.IndexType.*
import static cc.redberry.groovy.RedberryPhysics.DiracTrace
import static cc.redberry.groovy.RedberryPhysics.mandelstam
import static cc.redberry.groovy.RedberryStatic.*

/**
 * Created by poslavsky on 02/04/15.
 */
class SetupTest {
    @Before
    public void before() throws Exception {
        CC.reset()
    }

    @Test
    public void testInitialisation() {
        use(Redberry) {
            for (def cc in [true, false]) {
                CC.reset()//clear all previous definitions

                def defs = new Setup(cc)
                def t = 'Tr[V_aA*V_bB]'.t
                t <<= defs.FeynmanRules & defs.dTrace & defs.uTrace
                assert t == '-2*g**2*g_{ab}*g_{BA}'.t
            }
        }
    }

    @Test
    public void testTaylor() {
        use(Redberry) {
            def stp = new Setup(false)
            def t = stp.taylor('q_i') >> 'q_i*k^i'.t
            assert t == 'q_i*k^i'.t
        }
    }

    @Test
    public void testMathematica() throws Exception {
        use(Redberry) {
            def stp = new Setup(false, false)
            def t = '(Sin[x]*Sin[y] - Cos[x]*Cos[y])*g_mn + Cos[x+y]*g_mn'.t
            assert stp.mSimplify >> t == 0.t

            t = 'a*c + I*b*c - I*a*d + b*d'.t
            println stp.wFactor >> t
        }
    }

    @Test
    public void testMaple() throws Exception {

        use(Redberry) {
            def stp = new Setup(false, false)
            def t = 'x**2 - x'.t
            assert stp.mapleFactorTr >> t == 'x*(x - 1)'.t
        }
    }


    @Test
    public void testEffectiveVertex1() throws Exception {
        use(Redberry) {
            def stp = new Setup(true, true)
            def vertices = stp.effectiveQuarkoniaVertices()
            for (fl in ['charm', 'bottom'])
                for (spin in ['scalar', 'axial', 'tensor']) {
                    def t = vertices[spin] >> "A${spin}_{aA bB}[${fl}, k1_i, k2_i]".t
                    def expectedContent = [
                            'eps2^{b}[h2]', 'eps1^{a}[h1]',
                            "eps^{f}[h[$fl]]", "eps^{fg}[h[$fl]]", "m[$fl]", "p_{b}[$fl]",
                            'k1^{i}', 'k2^{i}',
                            'e_{abcd}', 'g_{BA}', 'g', 'g_{ab}']
                            .collect({ x -> x.t.hashCode() }).toSet()
                    def actualContent = TensorUtils.getAllDiffSimpleTensors(t).collect({ x -> x.t.hashCode() }).toSet()
                    assert expectedContent.containsAll(actualContent)
                }
        }
    }

    @Test
    public void testEffectiveVertex2() throws Exception {
        use(Redberry) {
            def stp = new Setup(false)
            def vertex = stp.effectivePairVertex()
            for (fl in ['fl', 'charm', 'bottom']) {
                def t = vertex >> "B_{aA bB}[${fl}, k1_i, k2_i]".t
                def expectedContent = [
                        'm[fl]', 'T_{A}', 'g', 'v[p2_{m}[fl]]',
                        'p1_{c}[fl]', 'cu[p1_{m}[fl]]',
                        'G_{b}', 'k1^{c}', 'k2_{m}', 'g_ab', 'g_AB'].t
                        .collect({ x -> x.t.hashCode() }).toSet()
                def actualContent = TensorUtils.getAllDiffSimpleTensors(t).collect({ x -> x.t.hashCode() }).toSet()
                assert expectedContent.containsAll(actualContent)
            }
        }
    }

    /**
     * Testing Ward identities
     */
    @Test
    public void testEffectiveVerticesWard() {
        use(Redberry) {

            def stp = new Setup(true)
            def vertices = stp.effectiveQuarkoniaVertices()

            for (def spin in ['scalar', 'axial', 'tensor']) {
                def vertex = vertices[spin]

                // Checking that A(g + g -> Q) is symmetric
                for (def f in ['charm', 'bottom'])
                    assert (vertex & ExpandAll) >> "A${spin}_{aA bB}[${f}, k1_m, k2_m] - A${spin}_{bB aA}[${f}, k2_m, k1_m]".t == 0.t

                // A(g + g -> Q) does not contains relative momentum
                def A = vertex >> "A${spin}_{aA bB}[fl, k1_m, k2_m]".t
                A.parentAfterChild { t ->
                    assert !('q_i[fl]'.t % t).exists
                }

                // Ward identities
                for (def temp in [
                        vertex >> "k1^a*A${spin}_{aA bB}[charm, k1_i, k2_i]".t,
                        vertex >> "k2^b*A${spin}_{aA bB}[charm, k1_i, k2_i]".t]) {
                    temp = ('k1_a = p_a[charm] - k2_a'.t
                            & stp.fullSimplify
                            & stp.massesSubs
                            & Factor
                    ) >> temp
                    assert temp == 0.t
                }
            }
        }
    }

    @Test
    public void testEffectivePairVertex1() {
        use(Redberry) {
            def stp = new Setup(false)
            def vertex = stp.effectivePairVertex()

            // Checking that A(g + g -> Q) is symmetric
            for (def f in ['charm', 'bottom'])
                assert (vertex & ExpandAll) >> "B_{aA bB}[${f}, k1_m, k2_m] - B_{bB aA}[${f}, k2_m, k1_m]".t == 0.t

        }
    }

    @Test
    public void testEffectivePairVertexWard() {
        use(Redberry) {
            def stp = new Setup(true)

            //two diagrams
            def g2 = 'cu[p1_m[fl]]*(V_aA*D[p1_m[fl] - k1_m, m[fl]]*V_bB + V_bB*D[p1_m[fl] - k2_m, m[fl]]*V_aA)*v[p2_m[fl]]'.t
            //third diagram (3-gluon)
            def g3 = 'cu[p1_m[fl]]*V^cC*v[p2_m[fl]]*G_c^d[k1_a + k2_a]*V_{aA bB dC}[k1_a, k2_a, -k1_a -k2_a]'.t
            def M = g2 + g3
            // Simplifying
            M <<= stp.FeynmanRules & ExpandAll[EliminateMetrics] & EliminateMetrics &
                    'p1_m[fl]*p1^m[fl] = m[fl]**2'.t & 'p2_m[fl]*p2^m[fl] = m[fl]**2'.t
            def pairVertex = 'B_{aA bB}[fl, k1_m, k2_m]'.t.eq(M)


            def mndlst = setMandelstam([k1_m: '0', k2_m: '0', 'p1_m[charm]': 'm[charm]', 'p2_m[charm]': 'm[charm]'])
            def simplifyMetrics = EliminateMetrics &
                    'eps1^a[h1] * k1_a = 0'.t &
                    'eps2^a[h2] * k2_a = 0'.t &
                    mndlst &
                    'd^i_i = 4'.t & 'd^A_A = 8'.t & "d^i'_i' = 4".t & "d^A'_A' = 3".t
            def fullSimplify = simplifyMetrics &
                    ExpandAll[simplifyMetrics] & simplifyMetrics &
                    stp.leviSimplify &
                    ExpandAll[simplifyMetrics] & simplifyMetrics

            def dTraceSimplify = DiracTrace[[Gamma: 'G_a', Simplifications: fullSimplify]]

            for (def amp in [
                    'k1^a * eps2^b[h2] * B_{aA bB}[charm, k1_m, k2_m]'.t,
                    'eps1^a[h1] * k2^b * B_{aA bB}[charm, k1_m, k2_m]'.t
            ]) {
                amp <<= pairVertex & fullSimplify

                def sum = new SumBuilder()
                for (int i = 0; i < amp.size(); ++i)
                    for (int j = 0; j < amp.size(); ++j) {
                        def ampC = amp[j]
                        ampC <<= Conjugate & Reverse[Matrix1, Matrix2] & stp.conjugateSpinors
                        ampC <<= { expr -> (expr.indices.free.si % expr.indices.free.si.inverted) >> expr } as Transformation

                        def amp2 = amp[i] * ampC
                        def indexless = amp2.indexlessSubProduct
                        def tensor = amp2.dataSubProduct

                        tensor <<= stp.epsSum & stp.uTrace & mndlst & dTraceSimplify &
                                fullSimplify & stp.massesSubs & stp.uSimplify
                        sum << indexless * tensor
                    }

                def amp2 = sum.build()
                amp2 <<= stp.massesSubs
                amp2 <<= 'u = 2*mc**2 - s - t'.t

                assert stp.wolframFactorTr >> amp2 == 0.t
            }
        }
    }

    /**
     * In case of axial meson, check Landau-Yang theorem
     */
    @Test
    public void testEffectiveVertexAxial() {
        use(Redberry) {
            def stp = new Setup(true)
            def vertices = stp.effectiveQuarkoniaVertices()
            def vertex = vertices['axial']
            assert ({
                def temp = vertex >> 'eps1^a[h1]*eps2^b[h2]*Aaxial_{aA bB}[bottom, k1_i, k2_i]'.t
                println TensorUtils.getAllDiffSimpleTensors(temp)
                ('p_a[bottom] = k1_a + k2_a'.t
                        & stp.fullSimplify
                        & 'm[bottom] = mb'.t.hold & 's = (2*mb)**2'.t
                        & Factor
                ) >> temp
            }()) == 0.t
        }
    }

    @Test
    public void testXYZ() throws Exception {
        use(Redberry) {
            def stp = new Setup(false, true)
            println stp.setXYZ('bottom', 'z', 'x')
        }
    }

    @Ignore
    @Test
    public void testMultiplicationTable() throws Exception {
        use(Redberry) {
            def stp = new SetupCC()
            def momentums = ['p_i[bottom]', 'p1_i[charm]', 'p2_i[charm]', 'k1_i', 'k2_i'].t

            def allExprs = []
            //single gamma
            momentums.each {
                def a = it
                allExprs << "g_AB*cu[p1_m[charm]]*G^i*$a*v[p2_m[charm]]".t
                allExprs << "cu[p1_m[charm]]*T_A*T_B*G^i*$a*v[p2_m[charm]]".t
                allExprs << "cu[p1_m[charm]]*T_B*T_A*G^i*$a*v[p2_m[charm]]".t
            }

            for (int i = 0; i < momentums.size(); ++i)
                for (int j = 0; j < momentums.size(); ++j) {
                    if (i == j)
                        continue
                    def a = momentums[i]
                    def b = '{_i -> _j}'.mapping >> momentums[j]
                    allExprs << "g_AB*cu[p1_m[charm]]*G^i*G^j*$a*$b*v[p2_m[charm]]".t
                    allExprs << "cu[p1_m[charm]]*T_A*T_B*G^i*G^j*$a*$b*v[p2_m[charm]]".t
                    allExprs << "cu[p1_m[charm]]*T_B*T_A*G^i*G^j*$a*$b*v[p2_m[charm]]".t
                }

            def tr = Identity

            println allExprs.size()
            def conjugate = Conjugate & InvertIndices
            conjugate &= Reverse[Matrix1, Matrix2]
            conjugate &= stp.conjugateSpinors

            //75
            //122876
            //DescriptiveStatistics:
            //n: 5625
            //min: 5.0
            //max: 144.0
            //mean: 21.844622222222245
            //std dev: 10.491243044688652
            //median: 21.0
            //skewness: 1.0270421570982076
            //kurtosis: 3.764432016727073
            //
            //
            //Process finished with exit code 0
            //661718
            //DescriptiveStatistics:
            //n: 5625
            //min: 74.0
            //max: 358.0
            //mean: 117.63875555555543
            //std dev: 33.90257764322778
            //median: 109.0
            //skewness: 2.2500673888050655
            //kurtosis: 7.4512299782792315

            DescriptiveStatistics raw = new DescriptiveStatistics()
            def totalRaw = 0
            for (int i = 0; i < allExprs.size(); ++i)
                for (int j = 0; j < allExprs.size(); ++j) {
                    def a = allExprs[i], b = allExprs[j]
                    b <<= conjugate
                    def tensor = a * b
                    def tm = timing({
                        tensor <<= stp.epsSum & stp.uTrace & stp.mandelstam & stp.dTraceSimplify &
                                stp.fullSimplify & stp.massesSubs & stp.uSimplify & stp.massesSubs & stp.wolframFactorTr
                    }, false)
                    raw.addValue(1.0 * tm)
                    totalRaw += tm

                    tr &= (a * b).eq tensor
                }
            println totalRaw
            println raw

            DescriptiveStatistics subs = new DescriptiveStatistics()
            def totalSubs = 0
            for (int i = 0; i < allExprs.size(); ++i)
                for (int j = 0; j < allExprs.size(); ++j) {
                    def a = allExprs[i], b = allExprs[j]
                    b <<= conjugate
                    def tensor = a * b
                    def tm = timing({
                        tensor <<= tr
                    }, false)
                    subs.addValue(1.0 * tm)
                    totalSubs += tm

                    if (!TensorUtils.isSymbolic(tensor))
                        println 'pizda'
                }

            println totalSubs
            println subs
        }

    }

    @Ignore
    @Test
    public void testMultiplicationTable1() throws Exception {
        use(Redberry) {
            def stp = new SetupCC()

            def subs = []
            //single gamma
            subs << "cu[p1_m[charm]]*T_A*T_B*v[p2_m[charm]]".t
            subs << "g_AB*cu[p1_m[charm]]*v[p2_m[charm]]".t
            subs << "g_AB*cu[p1_m[charm]]*G^i*v[p2_m[charm]]".t
            subs << "g_AB*cu[p1_m[charm]]*G^i*G5*v[p2_m[charm]]".t
            subs << "cu[p1_m[charm]]*T_A*T_B*G^i*v[p2_m[charm]]".t
            subs << "cu[p1_m[charm]]*T_A*T_B*G^i*G5*v[p2_m[charm]]".t
            subs << "g_AB*cu[p1_m[charm]]*G^i*G^j*v[p2_m[charm]]".t
            subs << "g_AB*cu[p1_m[charm]]*G^i*G^j*G5*v[p2_m[charm]]".t
            subs << "cu[p1_m[charm]]*T_A*T_B*G^i*G^j*v[p2_m[charm]]".t
            subs << "cu[p1_m[charm]]*T_A*T_B*G^i*G^j*G5*v[p2_m[charm]]".t

            subs << "f_ABC*cu[p1_m[charm]]*T^C*v[p2_m[charm]]".t
            subs << "d_ABC*cu[p1_m[charm]]*T^C*v[p2_m[charm]]".t


            for (int i = 0; i < subs.size(); ++i) {
                def l = "L${i + 1}${subs[i].indices.free}".t
                subs[i] = subs[i].eq l
            }

            def conjugate = Conjugate
            conjugate &= InvertIndices[LatinUpper]
            conjugate &= Reverse[Matrix1, Matrix2]
            conjugate &= stp.conjugateSpinors
            def ren = '{i -> a, j -> b}'.mapping

            for (int i = 0; i < subs.size(); ++i)
                for (int j = 0; j < subs.size(); ++j) {
                    def lhs = subs[i][1] * ((conjugate & ren) >> subs[j][1])

                    def rhs1 = subs[i][0] * ((conjugate & ren) >> subs[j][0])
                    rhs1 <<= stp.epsSum & stp.uTrace & stp.mandelstam & stp.dTraceSimplify &
                            stp.fullSimplify & stp.massesSubs & stp.uSimplify & stp.massesSubs & stp.wFactor

                    def rhs2 = (conjugate >> subs[i][0]) * (ren >> subs[j][0])
                    rhs2 <<= stp.epsSum & stp.uTrace & stp.mandelstam & stp.dTraceSimplify &
                            stp.fullSimplify & stp.massesSubs & stp.uSimplify & stp.massesSubs & stp.wFactor

                    println(stp.mFactor >> (rhs1 - rhs2))
                }
        }
    }

    @Test
    public void testName() throws Exception {
        use(Redberry) {
            def f = new File('/Users/poslavsky/Projects/redberry/redberry-pairedchi/output/cc1-1.redberry')
            StringBuilder sb = new StringBuilder()
            f.eachLine {
                sb.append(it)
            }
            def t = sb.toString().t

            t <<= 't1=-684062/1000'.t &
                    't2=-110326/1000'.t &
                    'u1=-434839/100000'.t &
                    'u2=-499841/10000'.t &
                    's=900'.t &
                    'mc = 15/10'.t &
                    'mb = 3'.t &
                    'g= 1'.t
            println t

            println TensorUtils.info(t)
        }
    }
}

