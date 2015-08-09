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
package cc.redberry.groovy.feyncalc.pairedchi.cli

import cc.redberry.core.context.CC
import cc.redberry.groovy.Redberry
import cc.redberry.groovy.feyncalc.pairedchi.SetupCC

import static cc.redberry.core.context.OutputFormat.SimpleRedberry
import static cc.redberry.core.context.OutputFormat.WolframMathematica
import static cc.redberry.core.utils.TensorUtils.info
import static cc.redberry.groovy.RedberryStatic.*

/**
 * Created by poslavsky on 28/07/15.
 */
class Main1 {
    public static void main(String[] args) {
        use(Redberry) {
            for (def S in ['x', 'y', 'z'])
                for (def L in ['x', 'y', 'z'])
                    for (def eps1 in [1, -1])
                        for (def eps2 in [1, -1]) {
                            CC.reset()
                            CC.resetTensorNames(-2907357143612290431)

                            SetupCC stp = new SetupCC()
                            //axial: 21369 (not projected), 51048 (1, 1)
                            def bottomSpin = 'scalar'

                            def pol = stp.setupPolarisations(eps1, eps2)
                            def qXYZ = stp.setXYZ('bottom', S, L)
                            def qpol = qXYZ['tr']


                            def diags = stp.diagrams(bottomSpin)


                            def i = 0
                            def file = new File("/Users/poslavsky/Downloads/amps_${eps1}_${eps2}_${S}_${L}.txt")
                            file.delete()
                            file << "Gluons factor (inversed):"
                            file << (stp.overallPolarizationFactor**(1 / 2)).toString(WolframMathematica)
                            file << "Bottom factor (inversed):"
                            file << qXYZ['den'].toString(WolframMathematica)


                            file << "\n\nMandelstam variables:\n\n"
                            stp.mandelstam.each {
                                file << "${it.toString(SimpleRedberry)};\n"
                            }

                            file << "\n\nSpinor structures:\n\n"
                            stp.spinorStructures.each {
                                file << "${it.toString(SimpleRedberry)};\n"
                            }

                            def tot = 0.t

                            def stringify = 'p1_a[charm] = p1_a'.t.hold &
                                    'p2_a[charm] = p2_a'.t.hold &
                                    'eps1_a[h1] = eps1_a'.t.hold &
                                    'eps2_a[h2] = eps2_a'.t.hold &
                                    'p_a[bottom] = p_a'.t.hold &
                                    'eps_{a}[h[bottom]] = epsP_a'.t.hold &
                                    'eps_{ab}[h[bottom]] = epsP_ab'.t.hold

                            file << "\n\nAmplitudes:\n\n"
                            for (def diag in diags) {
                                stp.log "amp $i of ${diags.size()}"
                                def amp = stp.calcAmplitude(diag, pol & qpol)
                                amp <<= stringify


                                def num = Numerator >> amp

//                def f = new File("/Users/poslavsky/Downloads/amp$i")
//                f.delete()
//                f << num.toString(OutputFormat.Redberry)

                                num <<= ExpandTensors & EliminateMetrics & stp.mFactor

                                amp = num / (Denominator >> amp)
                                tot += amp

                                amp <<= InvertIndices
                                file << "AMP[$i][A, B] = ${amp.toString(WolframMathematica)};\n\n"
                                ++i
                            }

                            println info(tot)
                        }
        }
    }
}
