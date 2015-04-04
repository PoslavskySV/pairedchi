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

import cc.redberry.core.tensor.Tensor
import cc.redberry.core.transformations.Transformation
import cc.redberry.core.utils.TensorUtils
import cc.redberry.groovy.Redberry

import static cc.redberry.core.context.OutputFormat.Maple
import static cc.redberry.core.context.OutputFormat.WolframMathematica
import static cc.redberry.groovy.RedberryStatic.*

/**
 * Created by poslavsky on 03/04/15.
 */
class ChiBCC extends Setup {
    private def qVertices, ccVertex

    private Map quarkDiagrams = [:]
    private Map gluonDiagrams = [:]
    private Tensor quarkDiagramsNotProjected = null;
    private Tensor gcc = null;

    private Transformation spinorsSimplify

    ChiBCC() {
        super(false, true);
        use(Redberry) {
            qVertices = effectiveQuarkoniaVertices().values() as Transformation
            ccVertex = effectivePairVertex()
            spinorsSimplify = Identity
            spinorsSimplify &= 'cu[p1_m[charm]]*G_a*p1^a[charm] = m[charm]*cu[p1_m[charm]]'.t
            spinorsSimplify &= 'G_a*p2^a[charm]*v[p2_m[charm]] = m[charm]*v[p2_m[charm]]'.t
        }
    }

    Tensor getGluonDiagrams(bottomSpin) {
        if (gluonDiagrams[bottomSpin] != null)
            return gluonDiagrams[bottomSpin]
        use(Redberry) {

            def simplify = FeynmanRules & qVertices & ccVertex & fullSimplify
            //first diagram
            def Ma = '1'.t
            Ma *= simplify >> "eps1^a[h1] * B_{aA cC}[charm, k1_i, -k1_i + p1_i[charm] + p2_i[charm]]".t
            Ma *= simplify >> 'G^cd[k1_i - p1_i[charm] - p2_i[charm]] * g^CD'.t
            Ma *= simplify >> "eps2^b[h2] * A${bottomSpin}_{dD bB}[bottom, -k2_i + p_i[bottom], k2_i]".t

            //second diagram
            def Mb = '1'.t
            Mb *= simplify >> "eps1^a[h1] * A${bottomSpin}_{aA cC}[bottom, k1_i, -k1_i + p_i[bottom]]".t
            Mb *= simplify >> 'G^cd[k1_i - p_i[bottom]] * g^CD'.t
            Mb *= simplify >> "eps2^b[h2] * B_{dD bB}[charm, -k2_i + p1_i[charm] + p2_i[charm], k2_i]".t

            def M = Ma + Mb
            M <<= fullSimplify & massesSubs & mFactor & spinorsSimplify & massesSubs & mFactor
            return (gluonDiagrams[bottomSpin] = M)
        }
    }


    Tensor getQuarkDiagramsNotProjected() {
        if (quarkDiagramsNotProjected != null)
            return quarkDiagramsNotProjected

        use(Redberry) {
            gcc = 'Vcc_iI = cu[p1_m[charm]]*G_i*T_I*v[p2_m[charm]]'.t

            quarkDiagramsNotProjected = '0'.t
            // (1,2,3)
            quarkDiagramsNotProjected += 'cu[p1_m[bottom]]*V_cC*Vcc^cC*D[p1_m[bottom] + pCharm_m, m[bottom]]*V_bB*eps2^b[h2]*D[k1_m - p2_m[bottom], m[bottom]]*V_aA*eps1^a[h1]*v[p2_m[bottom]]'.t
            // (3,2,1)
            quarkDiagramsNotProjected += 'cu[p1_m[bottom]]*V_aA*eps1^a[h1]*D[p1_m[bottom] - k1_m, m[bottom]]*V_bB*eps2^b[h2]*D[-pCharm_m - p2_m[bottom], m[bottom]]*V_cC*Vcc^cC*v[p2_m[bottom]]'.t
            // (1,3,2)
            quarkDiagramsNotProjected += 'cu[p1_m[bottom]]*V_bB*eps2^b[h2]*D[p1_m[bottom] - k2_m, m[bottom]]*V_cC*Vcc^cC*D[k1_m - p2_m[bottom], m[bottom]]*V_aA*eps1^a[h1]*v[p2_m[bottom]]'.t
            // (3,1,2)
            quarkDiagramsNotProjected += 'cu[p1_m[bottom]]*V_bB*eps2^b[h2]*D[p1_m[bottom] - k2_m, m[bottom]]*V_aA*eps1^a[h1]*D[-pCharm_m - p2_m[bottom], m[bottom]]*V_cC*Vcc^cC*v[p2_m[bottom]]'.t
            // (2,3,1)
            quarkDiagramsNotProjected += 'cu[p1_m[bottom]]*V_aA*eps1^a[h1]*D[p1_m[bottom] - k1_m, m[bottom]]*V_cC*Vcc^cC*D[k2_m - p2_m[bottom], m[bottom]]*V_bB*eps2^b[h2]*v[p2_m[bottom]]'.t
            // (2,1,3)
            quarkDiagramsNotProjected += 'cu[p1_m[bottom]]*V_cC*Vcc^cC*D[p1_m[bottom] + pCharm_m, m[bottom]]*V_aA*eps1^a[h1]*D[k2_m - p2_m[bottom], m[bottom]]*V_bB*eps2^b[h2]*v[p2_m[bottom]]'.t


            def masses = 'p2_{f}[bottom]*p2^{f}[bottom] = m[bottom]**2'.t & 'p1_{d}[bottom]*p1^{d}[bottom] = m[bottom]**2'.t
            quarkDiagramsNotProjected <<= 'pCharm_m = p1_m[charm] + p2_m[charm]'.t
            quarkDiagramsNotProjected <<= FeynmanRules & spinSingletProjector['bottom']
            quarkDiagramsNotProjected <<= dTraceSimplify & uTrace
            quarkDiagramsNotProjected <<= masses
            quarkDiagramsNotProjected <<= momentums['bottom']
            quarkDiagramsNotProjected <<= 'q_i[bottom] = q_i'.t.hold & taylor('q_i') & 'q_i = q_i[bottom]'.t
            quarkDiagramsNotProjected <<= mandelstam & ExpandAll[EliminateMetrics & mandelstam]
            return quarkDiagramsNotProjected
        }
    }

    Tensor getQuarkDiagrams(String bottomSpin) {
        if (quarkDiagrams[bottomSpin] != null)
            return quarkDiagrams[bottomSpin]
        use(Redberry) {
            log "Setting up quak diagrams for $bottomSpin ..."
            def Mc = getQuarkDiagramsNotProjected()
            Mc <<= totalSpinProjector[bottomSpin]
            Mc <<= fullSimplify & massesSubs & mFactor & gcc & spinorsSimplify & massesSubs & mFactor
            log "... done"
            return (quarkDiagrams[bottomSpin] = Mc)
        }
    }

    Tensor setupFeynmanDiagrams(bottomSpin) {
        use(Redberry) {


            def simplify = FeynmanRules & qVertices & ccVertex & fullSimplify
            //first diagram
            def Ma = '1'.t
            Ma *= simplify >> "eps1^a[h1] * B_{aA cC}[charm, k1_i, -k1_i + p1_i[charm] + p2_i[charm]]".t
            Ma *= simplify >> 'G^cd[k1_i - p1_i[charm] - p2_i[charm]] * g^CD'.t
            Ma *= simplify >> "eps2^b[h2] * A${bottomSpin}_{dD bB}[bottom, -k2_i + p_i[bottom], k2_i]".t

            //second diagram
            def Mb = '1'.t
            Mb *= simplify >> "eps1^a[h1] * A${bottomSpin}_{aA cC}[bottom, k1_i, -k1_i + p_i[bottom]]".t
            Mb *= simplify >> 'G^cd[k1_i - p_i[bottom]] * g^CD'.t
            Mb *= simplify >> "eps2^b[h2] * B_{dD bB}[charm, -k2_i + p1_i[charm] + p2_i[charm], k2_i]".t

            def Mc = getQuarkDiagrams(bottomSpin)

            def M = Ma + Mb + Mc
            M <<= 'cu[p1_a[charm]]*p1_a[charm]*G^a = -m[charm]*cu[p1_a[charm]]'.t &
                    'p2_a[charm]*G^a*v[p2_a[charm]] = m[charm]*v[p2_a[charm]]'.t

            return M
        }
    }

    Tensor calculateSquaredMatrixElement(String bottomSpin) {
        squareMatrixElement(setupFeynmanDiagrams(bottomSpin), "bottom: $bottomSpin")
    }

    void calculateAllProcesses(String output) {
        File file = new File(output)
        if (file.exists()) {
            file.delete();
            file = new File(output)
        }
        File fileMaple = new File(output)

        use(Redberry) {

            for (def bottomSpin in ['scalar', 'axial', 'tensor']) {
                def squared = calculateSquaredMatrixElement(bottomSpin)


                assert TensorUtils.isSymbolic(squared)

                def toStr = [scalar: 0, axial: 1, tensor: 2]
                def stringResult = (squared).toString(WolframMathematica)
                file << "chiB${toStr[bottomSpin]}cc = ${stringResult};"
                file << "\n"

                stringResult = (squared).toString(Maple)
                fileMaple << "chiB${toStr[bottomSpin]}cc := ${stringResult}:"
                fileMaple << "\n"
            }
        }
    }
}
