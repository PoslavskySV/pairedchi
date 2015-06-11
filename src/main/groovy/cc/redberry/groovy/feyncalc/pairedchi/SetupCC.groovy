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
import cc.redberry.core.tensor.Sum
import cc.redberry.core.tensor.Tensor
import cc.redberry.groovy.Redberry

import static cc.redberry.groovy.RedberryStatic.*

/**
 * Created by poslavsky on 03/04/15.
 */
class SetupCC extends Setup {
    SetupCC() {
        super(false, true);
    }

    def diagrams(bottomSpin) {
        use(Redberry) {
            def diagrams = []
            //gluon diagrams
            def glMa = """eps1^a[h1] * B_{aA cC}[charm, k1_i, k2_i - p_i[bottom]]
                        * G^cd[-k2_i + p_i[bottom]] * g^CD
                        * eps2^b[h2] * A${bottomSpin}_{dD bB}[bottom, -k2_i + p_i[bottom], k2_i]""".t
            def glMb = """eps2^b[h2] * B_{bB cC}[charm, k2_i, k1_i - p_i[bottom]]
                        * G^cd[-k1_i + p_i[bottom]] * g^CD
                        * eps1^a[h1] * A${bottomSpin}_{dD aA}[bottom, -k1_i + p_i[bottom], k1_i]""".t
            diagrams += [glMa, glMb]

            //3-gluon diagrams
            def gl3Ma = """eps1^a[h1] * Vcc^cC * V_{cC aA dD}[-p1_a[charm] - p2_a[charm], k1_a, k2_a - p_a[bottom]]
                        * G^de[-k2_a + p_a[bottom]]*g^DE
                        * A${bottomSpin}_{eE bB}[bottom, -k2_a + p_a[bottom], k2_a] * eps2^b[h2]""".t
            def gl3Mb = """eps2^b[h2] * Vcc^cC * V_{cC bB dD}[-p1_a[charm] - p2_a[charm], k2_a, k1_a - p_a[bottom]]
                        * G^ed[-k1_a + p_a[bottom]]*g^ED
                        * A${bottomSpin}_{eE aA}[bottom, -k1_a + p_a[bottom], k1_a] * eps1^a[h1]""".t
            //third diagram (s-chanel)
            def gl3Mc = """eps1^a[h1] * eps2^b[h2] * V_{aA bB cC}[k1_a, k2_a, -k1_a - k2_a]
                        * G^cd[k1_a + k2_a]*g^CD
                        * A${bottomSpin}_{dD eE}[bottom, k1_a + k2_a, -p1_a[charm] - p2_a[charm]] * Vcc^eE""".t
            diagrams += [gl3Ma, gl3Mb, gl3Mc]

            log 'Setting up quark diagrams ...'
            //quark diagrams
            // (1,2,3)
            def qMa = 'cu[p1_m[bottom]]*V_cC*Vcc^cC*D[p1_m[bottom] + pCharm_m, m[bottom]]*V_bB*eps2^b[h2]*D[k1_m - p2_m[bottom], m[bottom]]*V_aA*eps1^a[h1]*v[p2_m[bottom]]'.t
            // (3,2,1)
            def qMb = 'cu[p1_m[bottom]]*V_aA*eps1^a[h1]*D[p1_m[bottom] - k1_m, m[bottom]]*V_bB*eps2^b[h2]*D[-pCharm_m - p2_m[bottom], m[bottom]]*V_cC*Vcc^cC*v[p2_m[bottom]]'.t
            // (1,3,2)
            def qMc = 'cu[p1_m[bottom]]*V_bB*eps2^b[h2]*D[p1_m[bottom] - k2_m, m[bottom]]*V_cC*Vcc^cC*D[k1_m - p2_m[bottom], m[bottom]]*V_aA*eps1^a[h1]*v[p2_m[bottom]]'.t
            // (3,1,2)
            def qMd = 'cu[p1_m[bottom]]*V_bB*eps2^b[h2]*D[p1_m[bottom] - k2_m, m[bottom]]*V_aA*eps1^a[h1]*D[-pCharm_m - p2_m[bottom], m[bottom]]*V_cC*Vcc^cC*v[p2_m[bottom]]'.t
            // (2,3,1)
            def qMe = 'cu[p1_m[bottom]]*V_aA*eps1^a[h1]*D[p1_m[bottom] - k1_m, m[bottom]]*V_cC*Vcc^cC*D[k2_m - p2_m[bottom], m[bottom]]*V_bB*eps2^b[h2]*v[p2_m[bottom]]'.t
            // (2,1,3)
            def qMf = 'cu[p1_m[bottom]]*V_cC*Vcc^cC*D[p1_m[bottom] + pCharm_m, m[bottom]]*V_aA*eps1^a[h1]*D[k2_m - p2_m[bottom], m[bottom]]*V_bB*eps2^b[h2]*v[p2_m[bottom]]'.t

            def gcc = 'Vcc_iI = G_ij[p1_m[charm] + p2_m[charm]]*cu[p1_m[charm]]*V^j_I*v[p2_m[charm]]'.t
            for (def qM in [qMa, qMb, qMc, qMd, qMe, qMf]) {
                diagrams += calcQAmp(qM, bottomSpin)
            }
            diagrams = diagrams.collect({ it << (gcc & EliminateMetrics) })
            log '... done'
            return diagrams
        }
    }

    def calcQAmp(Tensor amp, bottomSpin) {
        use(Redberry) {
            def factor = Factor[[FactorScalars: true, FactorizationEngine: wolframFactorTr]]
            def mm = 'p2_{f}[bottom]*p2^{f}[bottom] = m[bottom]**2'.t & 'p1_{d}[bottom]*p1^{d}[bottom] = m[bottom]**2'.t
            amp <<= FeynmanRules & 'pCharm_m = p1_m[charm] + p2_m[charm]'.t & ExpandDenominator & EliminateMetrics & mm & mandelstam
            amp <<= spinSingletProjector['bottom'] & dTraceSimplify & mm
            amp <<= momentums['bottom'] & 'q_i[bottom] = q_i'.t.hold & taylor('q_i') & 'q_i = q_i[bottom]'.t
            amp <<= ExpandTensors[EliminateMetrics] & EliminateMetrics & mandelstam
            amp <<= totalSpinProjector[bottomSpin] & mandelstam
            amp <<= Together & factor

            return amp.class == Sum ? amp.toList() : [amp]
        }
    }

    public static void calc(bottomSpin, g1, g2, File output) {
        use(Redberry) {
            SetupCC stp = new SetupCC();
            def diags = stp.diagrams(bottomSpin)
            def pols = stp.setupPolarisations(g1, g2)
            def M2 = stp.calcProcess(diags, pols)
            output << M2.toString(OutputFormat.Redberry)
        }
    }
}
