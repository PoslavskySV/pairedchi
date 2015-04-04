package cc.redberry.groovy.feyncalc.pairedchi

import cc.redberry.groovy.Redberry
import org.junit.Test

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
        stp.calculateAll()
    }
}
