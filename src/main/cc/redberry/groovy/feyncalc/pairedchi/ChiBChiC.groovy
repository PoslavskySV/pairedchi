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

import static cc.redberry.core.context.OutputFormat.WolframMathematica

/**
 * Created by poslavsky on 02/04/15.
 */
class ChiBChiC extends Setup {
    private def qVertices

    ChiBChiC() {
        super(true, true);
        use(Redberry) {
            qVertices = effectiveQuarkoniaVertices().values() as Transformation
        }
    }

    Tensor setupFeynmanDiagrams(charmSpin, bottomSpin) {
        use(Redberry) {


            def simplify = FeynmanRules & qVertices & fullSimplify
            //first diagram
            def Ma = '1'.t
            Ma *= simplify >> "eps1^a[h1] * A${charmSpin}_{aA cC}[charm, k1_i, -k1_i + p_i[charm]]".t
            Ma *= simplify >> 'G^cd[k1_i - p_i[charm]] * g^CD'.t
            Ma *= simplify >> "eps2^b[h2] * A${bottomSpin}_{dD bB}[bottom, -k2_i + p_i[bottom], k2_i]".t

            //second diagram
            def Mb = '1'.t
            Mb *= simplify >> "eps1^a[h1] * A${bottomSpin}_{aA cC}[bottom, k1_i, -k1_i + p_i[bottom]]".t
            Mb *= simplify >> 'G^cd[k1_i - p_i[bottom]] * g^CD'.t
            Mb *= simplify >> "eps2^b[h2] * A${charmSpin}_{dD bB}[charm, -k2_i + p_i[charm], k2_i]".t

            //total matrix element
            Ma + Mb
        }
    }

    void calculateAll() {

        File file = new File('/Users/poslavsky/Projects/redberry/redberry-groovy-scripts/src/main/groovy/cc/redberry/groovy/scripts/feyncalc/qcd/pairedChi/ChiCChiBResult.m')
        if (file.exists()) {
            file.delete();
            file = new File('/Users/poslavsky/Projects/redberry/redberry-groovy-scripts/src/main/groovy/cc/redberry/groovy/scripts/feyncalc/qcd/pairedChi/ChiCChiBResult.m')
        }

        use(Redberry) {
            for (def charmSpin in ['scalar', 'axial', 'tensor'])
                for (def bottomSpin in ['scalar', 'axial', 'tensor']) {
                    def diagram = setupFeynmanDiagrams(charmSpin, bottomSpin)
                    def squared = squareMatrixElement(diagram, "charm: $charmSpin, bottom: $bottomSpin")

                    //total spin projection
                    'sqrt[x] := x**(1/2)'.t
                    'gf#scalar[fl] := 1/sqrt[3]'.t
                    'gf#axial[fl] := 1/2/sqrt[2]/m[fl]'.t
                    'gf#tensor[fl] := 1'.t
                    //wave function and color
                    def wCoeff = 'R[fl] * sqrt[3/4/pi] * 2/sqrt[2]/sqrt[2*m[fl]] * sqrt[1/3]'.t
                    //spin projection
                    def sCoeff = '1/2/sqrt[2]/m[fl]'.t
                    //total coefficient
                    "gf[fl] := $wCoeff * $sCoeff".t
                    //cross section
                    'crs := 1/(16*pi*s**2)/(8*8*2*2)'.t

                    squared *= "crs * (gf[charm] * gf[bottom] * gf#${charmSpin}[charm] * gf#${bottomSpin}[bottom])**2".t
                    squared <<= massesSubs

                    assert TensorUtils.isSymbolic(squared)

                    def toStr = [scalar: 0, axial: 1, tensor: 2]
                    def stringResult = (wolframFactorTr >> squared).toString(WolframMathematica)
                    file << "chiC${toStr[charmSpin]}chiB${toStr[bottomSpin]} = ${stringResult};"
                    file << "\n"
                }
        }
    }
}
