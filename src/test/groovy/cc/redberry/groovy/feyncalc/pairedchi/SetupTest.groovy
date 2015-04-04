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
import cc.redberry.core.utils.TensorUtils
import cc.redberry.groovy.Redberry
import org.junit.Before
import org.junit.Test

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
                t <<= defs.dTrace & defs.uTrace

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
    public void testEffectiveVertex1() throws Exception {
        use(Redberry) {
            def stp = new Setup(true)
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
    public void testEffectiveVertex3() {
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

    /**
     * Testing Ward identities
     */
    @Test
    public void testEffectiveVertex4() {
        use(Redberry) {

            def stp = new Setup(false)
            def vertex = stp.effectivePairVertex()

            // Checking that A(g + g -> Q) is symmetric
            for (def f in ['charm', 'bottom'])
                assert (vertex & ExpandAll) >> "B_{aA bB}[${f}, k1_m, k2_m] - B_{bB aA}[${f}, k2_m, k1_m]".t == 0.t

            // Ward identities
            def tr = stp.fullSimplify &
                    'cu[p1_a[charm]]*p1_a[charm]*G^a = m[charm]*cu[p1_a[charm]]'.t &
                    'p2_a[charm]*G^a*v[p2_a[charm]] = -m[charm]*v[p2_a[charm]]'.t &
                    'v[p2_a[charm]] * cu[p1_a[charm]] = 1'.t &
                    stp.dTrace & stp.uTrace & stp.fullSimplify &
                    stp.massesSubs &
                    Factor

            def temp

            temp = vertex >> "k1^a*B_{aA bB}[charm, k1_i, k2_i]".t
            temp <<= 'k1_a = p1_a[charm] + p2_a[charm] - k2_a'.t & tr
            assert temp == 0.t

            temp = vertex >> "k2^b*B_{aA bB}[charm, k1_i, k2_i]".t
            temp <<= 'k2_a = p1_a[charm] + p2_a[charm] - k1_a'.t & tr
            assert temp == 0.t
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
    public void testPerformance() throws Exception {
        use(Redberry) {
            def stp = new Setup(false, true)
            def exprs = []
            Setup.class.classLoader.getResourceAsStream('expressions').eachLine { l ->
                exprs << l.t
            }

            println exprs[0]
            timing {
                for (int i = 0; i < exprs.size(); i++) {
                    exprs[i] <<= stp.epsSum & stp.uTrace & stp.dTraceSimplify &
                            stp.fullSimplify & stp.massesSubs
                    assert TensorUtils.isSymbolic(exprs[i])
                }
            }

        }
    }
}
