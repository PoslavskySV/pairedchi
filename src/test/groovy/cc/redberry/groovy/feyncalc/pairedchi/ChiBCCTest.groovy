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
import cc.redberry.core.tensor.Product
import cc.redberry.core.utils.TensorUtils
import cc.redberry.groovy.Redberry
import org.junit.Test

import static cc.redberry.core.utils.TensorUtils.getAllDiffSimpleTensors
import static cc.redberry.groovy.RedberryStatic.EliminateMetrics
import static cc.redberry.groovy.RedberryStatic.ExpandAll

/**
 * Created by poslavsky on 03/04/15.
 */
class ChiBCCTest {

    @Test
    public void testWard1() throws Exception {
        use(Redberry) {
            ChiBCC stp = new ChiBCC();

        }
    }

    @Test
    public void testBottomLeg() throws Exception {
        use(Redberry) {
            ChiBCC stp = new ChiBCC();

            def Mc = '0'.t
            def gcc = 'Vcc_iI = cu[p1_m[charm]]*G_i*T_I*v[p2_m[charm]]'.t
            // (1,2,3)
            Mc += 'cu[p1_m[bottom]]*V_cC*Vcc^cC*D[p1_m[bottom] + pCharm_m, m[bottom]]*V_bB*eps2^b[h2]*D[k1_m - p2_m[bottom], m[bottom]]*V_aA*eps1^a[h1]*v[p2_m[bottom]]'.t
            // (3,2,1)
            Mc += 'cu[p1_m[bottom]]*V_aA*eps1^a[h1]*D[p1_m[bottom] - k1_m, m[bottom]]*V_bB*eps2^b[h2]*D[-pCharm_m - p2_m[bottom], m[bottom]]*V_cC*Vcc^cC*v[p2_m[bottom]]'.t
            // (1,3,2)
            Mc += 'cu[p1_m[bottom]]*V_bB*eps2^b[h2]*D[p1_m[bottom] - k2_m, m[bottom]]*V_cC*Vcc^cC*D[k1_m - p2_m[bottom], m[bottom]]*V_aA*eps1^a[h1]*v[p2_m[bottom]]'.t
            // (3,1,2)
            Mc += 'cu[p1_m[bottom]]*V_bB*eps2^b[h2]*D[p1_m[bottom] - k2_m, m[bottom]]*V_aA*eps1^a[h1]*D[-pCharm_m - p2_m[bottom], m[bottom]]*V_cC*Vcc^cC*v[p2_m[bottom]]'.t
            // (2,3,1)
            Mc += 'cu[p1_m[bottom]]*V_aA*eps1^a[h1]*D[p1_m[bottom] - k1_m, m[bottom]]*V_cC*Vcc^cC*D[k2_m - p2_m[bottom], m[bottom]]*V_bB*eps2^b[h2]*v[p2_m[bottom]]'.t
            // (2,1,3)
            Mc += 'cu[p1_m[bottom]]*V_cC*Vcc^cC*D[p1_m[bottom] + pCharm_m, m[bottom]]*V_aA*eps1^a[h1]*D[k2_m - p2_m[bottom], m[bottom]]*V_bB*eps2^b[h2]*v[p2_m[bottom]]'.t

            Mc <<= 'pCharm_m = p1_m[charm] + p2_m[charm]'.t

            def masses = 'p2_{f}[bottom]*p2^{f}[bottom] = m[bottom]**2'.t & 'p1_{d}[bottom]*p1^{d}[bottom] = m[bottom]**2'.t
            println 'xxxxxx'
            Mc <<= stp.FeynmanRules & stp.spinSingletProjector['bottom']
            Mc <<= stp.dTraceSimplify & stp.uTrace
            Mc <<= masses


            println Mc.size()
            println Mc[0]
            println Mc[1]
            println Mc[2]
            println getAllDiffSimpleTensors(Mc)

            Mc <<= stp.momentums['bottom']
            // Taylor expansion (scalar meson)
            Mc <<= 'q_i[bottom] = q_i'.t.hold & stp.taylor('q_i') & 'q_i = q_i[bottom]'.t

            println Mc.size()
            Mc <<= stp.mandelstam & ExpandAll[EliminateMetrics & stp.mandelstam] &
                    stp.totalSpinProjector['tensor']

            println Mc.size()
            Mc <<= stp.fullSimplify & stp.massesSubs & stp.mFactor & gcc

            println Mc.size()
            println Mc[0]
            println Mc[1]
        }

    }

    @Test
    public void test1123e4() throws Exception {
        use(Redberry) {
            ChiBCC stp = new ChiBCC();
            def t = stp.getGluonDiagrams('scalar')
            println t[0]
        }
    }

    @Test
    public void test1123() throws Exception {
        use(Redberry) {
            ChiBCC stp = new ChiBCC();
            for (def spin in ['scalar', 'axial', 'tensor']) {

                println 'GLUON'
                def t = stp.getGluonDiagrams(spin)
                println t.size()
                println t[0]
                println t[1]
                println t[2]

                println 'QUARK'
                t = stp.getQuarkDiagrams(spin)
                println t.size()
                println t[0]
                println t[1]
                println t[2]
                //todo assert content!
            }
        }
    }

    @Test
    public void test12323() throws Exception {
        use(Redberry) {

            def file = new File('/Users/poslavsky/Projects/redberry/redberry-groovy-scripts/src/main/groovy/cc/redberry/groovy/scripts/feyncalc/qcd/pairedChi/ChiBCCResult2.m')
            if (file.exists()) {
                file.delete()
                file = new File('/Users/poslavsky/Projects/redberry/redberry-groovy-scripts/src/main/groovy/cc/redberry/groovy/scripts/feyncalc/qcd/pairedChi/ChiBCCResult2.m')
            }
            ChiBCC stp = new ChiBCC();
            //for (def spin in ['scalar', 'axial', 'tensor']) {

//            def t =
            def r1 = stp.squareMatrixElement(stp.getQuarkDiagrams('tensor'))
            def r2 = stp.squareMatrixElement(stp.getGluonDiagrams('tensor'))
            def r = r1 + r2

            println TensorUtils.hasImaginaryPart(r)
            println r.size()
            println 'r1'
            //r1 <<= stp.wolframFactorTr
            //println r1

            println 'r2'
            //r2 <<= stp.wolframFactorTr
            //println r2


            file << ('r1=' + r1.toString(OutputFormat.WolframMathematica) + ';')
            file << '\n'

            file << ('r2=' + r2.toString(OutputFormat.WolframMathematica) + ';')
            file << '\n'
        }
    }

    @Test
    public void test123() throws Exception {
        use(Redberry) {
            Setup s = new Setup(false)
            s.mandelstam.transformations.each { expr ->
                Product lhs = expr[0]
                def rhs = s.massesSubs >> expr[1]
                def f = { x -> "Momentum[${x.stringName}]" }
                def str = ''
                str += "Pair[${f(lhs[0])},${f(lhs[1])}]"
                str += '='
                str += rhs.toString(OutputFormat.WolframMathematica)
                str += ';'
                println str
            }
        }
    }

    @Test
    public void testWardIdentities() throws Exception {
        use(Redberry) {
            ChiBCC stp = new ChiBCC()
            for (def charmSpin in ['scalar', 'axial', 'tensor'])
                for (def bottomSpin in ['scalar', 'axial', 'tensor']) {
                    def M = stp.setupFeynmanDiagrams(charmSpin, bottomSpin)
                    for (def i in [1, 2])
                        assert ({
                            def ward = ("eps${i}_i[h${i}] = k${i}_i".t
                                    & "k${3 - i}_i = p_i[bottom] - k${i}_i + p_i[charm]".t
                            ) >> M
                            ward <<= stp.fullSimplify & stp.massesSubs & stp.mFactor
                            ward
                        }()) == 0.t
                }
        }
    }

    @Test
    public void testCalculateAll() throws Exception {
        def stp = new ChiBCC()
        stp.calculateAll()
    }

    @Test
    public void test2() throws Exception {
        use(Redberry) {
            def stp = new ChiBCC()
            def mc = stp.mc('axial')

            println mc
        }
    }
}
