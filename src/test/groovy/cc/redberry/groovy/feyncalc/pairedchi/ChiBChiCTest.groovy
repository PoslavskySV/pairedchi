package cc.redberry.groovy.feyncalc.pairedchi

import cc.redberry.core.utils.TensorUtils
import cc.redberry.groovy.Redberry
import org.junit.Test

import static cc.redberry.core.context.OutputFormat.WolframMathematica

/**
 * Created by poslavsky on 02/04/15.
 */
class ChiBChiCTest {
    @Test
    public void testWardIdentities() throws Exception {
        use(Redberry) {
            ChiBChiC stp = new ChiBChiC()
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
        def stp = new ChiBChiC()

        File file = new File('/Users/poslavsky/Projects/redberry/redberry-groovy-scripts/src/main/groovy/cc/redberry/groovy/scripts/feyncalc/qcd/pairedChi/ChiCChiBResult.m')
        if (file.exists()) {
            file.delete();
            file = new File('/Users/poslavsky/Projects/redberry/redberry-groovy-scripts/src/main/groovy/cc/redberry/groovy/scripts/feyncalc/qcd/pairedChi/ChiCChiBResult.m')
        }

        use(Redberry) {
            for (def charmSpin in ['scalar', 'axial', 'tensor'])
                for (def bottomSpin in ['scalar', 'axial', 'tensor']) {
                    def diagram = stp.setupFeynmanDiagrams(charmSpin, bottomSpin)
                    def squared = stp.squareMatrixElement(diagram, "charm: $charmSpin, bottom: $bottomSpin")

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
                    squared <<= stp.massesSubs

                    assert TensorUtils.isSymbolic(squared)

                    def toStr = [scalar: 0, axial: 1, tensor: 2]
                    def stringResult = (stp.wolframFactorTr >> squared).toString(WolframMathematica)
                    file << "chiC${toStr[charmSpin]}chiB${toStr[bottomSpin]} = ${stringResult};"
                    file << "\n"
                }
        }
    }


    @Test
    public void test1() throws Exception {
        use(Redberry) {

        }
    }
}
