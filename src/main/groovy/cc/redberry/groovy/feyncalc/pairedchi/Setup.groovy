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
import cc.redberry.core.number.Complex
import cc.redberry.core.tensor.*
import cc.redberry.core.transformations.Transformation
import cc.redberry.core.utils.IteratorWithProgress
import cc.redberry.groovy.Redberry
import com.maplesoft.openmaple.Engine
import com.maplesoft.openmaple.EngineCallBacksDefault
import com.wolfram.jlink.KernelLink
import com.wolfram.jlink.MathLinkFactory

import static cc.redberry.core.context.OutputFormat.Maple
import static cc.redberry.core.context.OutputFormat.WolframMathematica
import static cc.redberry.core.indices.IndexType.Matrix1
import static cc.redberry.core.indices.IndexType.Matrix2
import static cc.redberry.core.utils.TensorUtils.info
import static cc.redberry.core.utils.TensorUtils.isSymbolic
import static cc.redberry.groovy.RedberryPhysics.*
import static cc.redberry.groovy.RedberryStatic.*

/**
 * Created by poslavsky on 02/04/15.
 */
class Setup implements AutoCloseable {
    public static final Map spins = [scalar: 0, axial: 1, tensor: 2]

    /**
     * Project c cbar onto charmonium or not
     */
    final boolean projectCC
    final long startTime
    final boolean log
    final boolean projectXYZ

    /** Feynman rules and auxiliary tensors */
    public def J2, J4, PS, V, D, G, V3, FeynmanRules

    /** Spin projectors and polarisation sums */

    /**
     *  Map flavour -> substitution
     */
    public Map spinSingletProjector
    /**
     * Map total spin -> substitution
     */
    public Map totalSpinProjector
    /**
     * Sum over polarisations
     */
    public Transformation epsSum

    /** Transformations */
    public def mandelstam, massesSubs, momentums,
               dTrace, dTraceSimplify, diracSimplify, spinorsSimplify, uTrace, uSimplify, leviSimplify,
               simplifyMetrics, fullSimplify, fullSimplifyE,
               conjugateSpinors, momentumConservation

    /**
     * Polarisations
     */
    public Tensor overallPolarizationFactor

    //Mathematica related
    KernelLink mathematicaKernel
    public Transformation wFactor, wolframFactorTr, mSimplify, wolframSimplifyTr, wolframFactorSqrtTr

    //Maple related
    //only one instance per programm (maple is slag)
    static volatile Engine mapleEngine;
    public Transformation mFactor, mapleFactorTr

    public Setup(boolean projectCC) {
        this(projectCC, false, false)
    }

    public Setup(boolean projectCC, boolean projectXYZ) {
        this(projectCC, projectXYZ, false)
    }

    public Setup(boolean projectCC, boolean projectXYZ, boolean log) {
        this.projectXYZ = projectXYZ
        this.startTime = System.currentTimeMillis()
        this.log = log
        this.projectCC = projectCC
        init()
    }

    void init() {//basic initialization
        use(Redberry) {
            log "Started with seed ${CC.nameManager.seed}"
            log 'Setting up basic tensors'

            //Setting up matrix objects
            defineMatrices 'T_A', Matrix2.matrix,//unitary matrices
                    'G_a', 'G5', Matrix1.matrix, //gamma matrices
                    'v[p_a]', 'u[p_a]', Matrix1.vector, Matrix2.vector, //final quark wave function
                    'cu[p_a]', 'cv[p_a]', Matrix1.covector, Matrix2.covector, //final antiquark wave function
                    'V_iA', Matrix1.matrix, Matrix2.matrix, //quark-gluon vertex
                    'D[p_m, mass]', Matrix1.matrix //quark propagator

            //Levi-Civita, SU(N) symmetric and structure tensors
            Quiet { Tensors.setAntiSymmetric('e_abcd') }
            Quiet { Tensors.setSymmetric('d_ABC') }
            Quiet { Tensors.setAntiSymmetric('f_ABC') }
            //polarization tensor
            Quiet { Tensors.setSymmetric('eps_ab[h]') }

            setupFeynRules()
            setupProjectors()
            setupTransformations()
            setupMathematica()
            setupPolarizationFactor()
            if (!projectCC)
                setupSpinorStructures()
            log 'Setup finished'
        }
    }

    Tensor wolframFunc(String command, Map bindings) {
        use(Redberry) {
            bindings.each { k, v ->
                command = command.replace(k, (v.class == Expression ? v[1] : v).toString(WolframMathematica))
            }
            mathematicaKernel.evaluateToInputForm(command, 0).replace('^', '**').t
        }
    }

    OutputStream dummyOut = new PrintStream(new OutputStream() {
        @Override
        void write(int b) throws IOException {
        }
    });

    Tensor mapleFunc(String command, Map bindings) {
        use(Redberry) {
            bindings.each { k, v ->
                def rhs = (v.class == Expression ? v[1] : v)
                command = command.replace(k, rhs.toString(Maple))
            }
            def out = System.out
            System.out = dummyOut;
            try {
                return mapleEngine.evaluate(command).toString().replace('^', '**').t
            } catch (Exception e) {
                System.out = out;
                throw new RuntimeException(e)
            } finally {
                System.out = out;
            }
        }
    }

    void setupMathematica() {
        log 'Setting up Mathematica and Maple'
        String[] args;
        switch (os()) {
            case 'mac':
                args = ["-linkmode", "launch", "-linkname", "\"/Applications/Mathematica.app/Contents/MacOS/MathKernel\" -mathlink"];
                break
            case 'linux':
                args = ["-linkmode", "launch", "-linkname", "math -mathlink"]
                break
            default:
                throw new RuntimeException('No Windows')
        }

        mathematicaKernel = MathLinkFactory.createKernelLink(args)
        mathematicaKernel.discardAnswer();
        if (mapleEngine == null)
            mapleEngine = new Engine(['java'] as String[], new EngineCallBacksDefault(), null, null)

        def mop = { tensor, func, bindings, e = this.&wolframFunc ->
            def t = e(func, bindings)
            return tensor.class == Expression ? tensor[0].eq(t) : t
        }
        wolframFactorTr = { tensor -> mop(tensor, 'Factor[expr]', [expr: tensor]) } as Transformation
        wolframSimplifyTr = { tensor -> mop(tensor, 'Simplify[expr, Reals]', [expr: tensor]) } as Transformation
        wolframFactorSqrtTr = { tensor -> mop(tensor, '(expr)/.{Power[x_, y_] :> Power[Factor[Expand[x]], y]}', [expr: tensor]) } as Transformation
        mapleFactorTr = { tensor -> mop(tensor, 'factor(expr);', [expr: tensor], this.&mapleFunc) } as Transformation

        wFactor = Factor[[FactorScalars: false, FactorizationEngine: wolframFactorTr]]
        mFactor = Factor[[FactorScalars: false, FactorizationEngine: mapleFactorTr]]

        mSimplify = Factor[[FactorScalars: false, FactorizationEngine: wolframSimplifyTr]]
    }

    /**
     * Feynman rules
     */
    void setupFeynRules() {
        log 'Setting up Feynman rules'
        use(Redberry) {
            PS = 'PS_mn[k_a] = -g_mn'.t
//            if (projectCC)
//                PS = 'PS_mn[k_a] = -g_mn'.t
//            else {
//                PS = 'PS_mn[k_a] = -g_mn - (k_m*n_n + k_n*n_m)/(k_a*n^a) + n_a*n^a*k_m*k_n/(k_a*n^a)**2'.t
//                //particular n
//                PS <<= 'n_a = p_a[bottom]'.t
//            }
            //Quark vertex
            V = 'V_mA = I*g*G_m*T_A'.t
            //Quark propagator
            D = 'D[p_m, mass] = I*(mass + p_m*G^m)/(p_m*p^m - mass**2)'.t
            //Gluon propagator
            G = PS >> 'G_mn[k_a] = I*PS_mn[k_a]/(k_a*k^a)'.t
            //3-gluon vertex
            V3 = 'V_{Aa Bb Cc}[k1_m, k2_m, k3_m] = g*f_{ABC}*(g_ab*(k1_c - k2_c) + g_bc*(k2_a - k3_a) + g_ca*(k3_b - k1_b))'.t
            //All Feynman rules in one transformation
            FeynmanRules = V & D & G & V3
        }
    }

    /**
     * Spin and polarisation sums
     */
    void setupProjectors() {
        log 'Setting up spin projectors and polarisation sums'
        use(Redberry) {
            // Auxiliary tensors
            J2 = 'J_ab[p_i, m] := p_a*p_b/m**2 - g_ab'.t
            J4 = 'J_abcd[p_i, m] := (J_ac[p_i, m]*J_bd[p_i, m] + J_ad[p_i, m]*J_bc[p_i, m])/2 - J_ab[p_i, m]*J_cd[p_i, m]/3'.t
            J4 <<= Expand

            /** quarkonia polarization projector */
            totalSpinProjector = [scalar: Identity, axial: Identity, tensor: Identity]
            /** sum over polarizations */

            epsSum = Identity
            conjugateSpinors = Identity
            spinSingletProjector = [:]
            momentums = [:]
            for (def fl in ['fl', 'charm', 'bottom']) { //for each flavour (`fl` denotes abstract flavour)
                //quarks momentums in terms of quarkonia momentum (q_i is relative momentum)
                momentums[fl] = "p2_m[$fl] = p_m[$fl]/2 + q_m[$fl]".t.hold & "p1_m[$fl] = p_m[$fl]/2 - q_m[$fl]".t.hold

                // set spin singlet projection
                spinSingletProjector[fl] = "v[p2_m[$fl]]*cu[p1_m[$fl]] = (p2_m[$fl]*G^m - m[$fl]) * epsS_m[$fl] * G^m * (p1_m[$fl]*G^m + m[$fl])".t

                // total spin projection
                if (projectXYZ) {
                    def s = "q_i[$fl]*epsS_j[$fl] = eps_i[hS[$fl]]*eps_j[hL[$fl]]".t
                    totalSpinProjector['scalar'] &= s
                    totalSpinProjector['axial'] &= s
                    totalSpinProjector['tensor'] &= s
                } else {
                    totalSpinProjector['scalar'] &= "q_i[$fl]*epsS_j[$fl] = -J_ij[p_i[$fl], 2*m[$fl]]".t
                    totalSpinProjector['axial'] &= "q_i[$fl]*epsS_j[$fl] = e_ijab * p^a[$fl] * eps^b[h[$fl]]".t
                    totalSpinProjector['tensor'] &= "q_i[$fl]*epsS_j[$fl] = eps_ij[h[$fl]]".t
                }

                //sum over polarizations for axial meson
                epsSum &= "eps_a[h[$fl]]*eps_b[h[$fl]] = J_ab[p_a[$fl], 2*m[$fl]]".t
                //sum over polarizations for tensor meson
                epsSum &= "eps_ab[h[$fl]]*eps_cd[h[$fl]] = J_abcd[p_a[$fl], 2*m[$fl]]".t

                // spinor sums
                epsSum &= "u[p1_m[$fl]]*cu[p1_m[$fl]] = m[$fl] + p1^m[$fl]*G_m".t
                epsSum &= "v[p2_m[$fl]]*cv[p2_m[$fl]] = -m[$fl] + p2^m[$fl]*G_m".t

                conjugateSpinors &= "v[p2_m[$fl]]*cu[p1_m[$fl]] = u[p1_m[$fl]]*cv[p2_m[$fl]]".t
            }

            //sum over polarizations for gluons
            epsSum &= PS >> 'eps1_a[h1]*eps1_b[h1] = PS_ab[k1_a]'.t
            epsSum &= PS >> 'eps2_a[h2]*eps2_b[h2] = PS_ab[k2_a]'.t
        }
    }

    /**
     * Setup all transformations
     */
    void setupTransformations() {
        log 'Setting up transformations'
        use(Redberry) {
            //Dirac, unitary traces and relative simplifications
            dTrace = DiracTrace
            uTrace = UnitaryTrace['T_A', 'f_ABC', 'd_ABC', '3']
            uSimplify = UnitarySimplify['T_A', 'f_ABC', 'd_ABC', '3']
            leviSimplify = LeviCivitaSimplify.minkowski['e_abcd']

            //Mandelstam variables and mass shell
            if (projectCC)
                mandelstam = setMandelstam([k1_m: '0', k2_m: '0', 'p_m[charm]': '2*m[charm]', 'p_m[bottom]': '2*m[bottom]'])
            else {
                mandelstam = setMandelstam5([k1_m: '0', k2_m: '0', 'p1_m[charm]': 'm[charm]', 'p2_m[charm]': 'm[charm]', 'p_m[bottom]': '2*m[bottom]'])
            }

            /** Simplifications of polarizations */
            def simplifyPolarizations = 'eps1^a[h1] * k1_a = 0'.t &
                    'eps2^a[h2] * k2_a = 0'.t &
                    'eps^a[h[charm]] * p_a[charm] = 0'.t &
                    'eps^a[h[bottom]] * p_a[bottom] = 0'.t &
                    'eps^a[hS[bottom]] * p_a[bottom] = 0'.t &
                    'eps^a[hL[bottom]] * p_a[bottom] = 0'.t &
                    'eps^ab[h[charm]] * p_a[charm] = 0'.t &
                    'eps^ab[h[bottom]] * p_a[bottom] = 0'.t

            //substituting masses for factor
            massesSubs = 'm[charm] = mc'.t.hold & 'm[bottom] = mb'.t.hold

            /** Full simplification */
            simplifyMetrics = EliminateMetrics & simplifyPolarizations & mandelstam &
                    'd^i_i = 4'.t & 'd^A_A = 8'.t & "d^i'_i' = 4".t & "d^A'_A' = 3".t
            if (projectCC)
                simplifyMetrics &= 'e_abcd * p^a[charm] * p^b[bottom] * k1^c * k2^d = 0'.t
            simplifyMetrics &= uSimplify

            fullSimplify = simplifyMetrics &
                    ExpandAll[simplifyMetrics] & simplifyMetrics &
                    leviSimplify &
                    ExpandAll[simplifyMetrics] & simplifyMetrics

            fullSimplifyE = simplifyMetrics &
                    ExpandTensors[simplifyMetrics] & simplifyMetrics &
                    leviSimplify &
                    ExpandTensors[simplifyMetrics] & simplifyMetrics

            dTraceSimplify = DiracTrace[[Simplifications: fullSimplify]]
            diracSimplify = DiracSimplify[[Simplifications: simplifyMetrics]]
            diracSimplify &= 'G_a*G_b*eps^ab[h[bottom]] = 0'.t
            diracSimplify &= 'G_a*G5*G_b*eps^ab[h[bottom]] = 0'.t
            diracSimplify &= 'G_a*G_c*G_b*eps^ab[h[bottom]] = 2*G_b*eps_c^b[h[bottom]]'.t

            spinorsSimplify = Identity
            spinorsSimplify &= SpinorsSimplify[[uBar: 'cu[p1_m[charm]]'.t, momentum: 'p1_m[charm]', mass: 'mc']]
            spinorsSimplify &= SpinorsSimplify[[v: 'v[p2_m[charm]]'.t, momentum: 'p2_m[charm]', mass: 'mc']]

            if (projectCC)
                momentumConservation = 'p_i[bottom] = k1_i + k2_i - p_i[charm]'.t.hold
            else
                momentumConservation = 'p_i[bottom] = k1_i + k2_i - p1_i[charm] - p2_i[charm]'.t.hold
        }
    }

    /**
     * First term in Taylor series
     * @param momentum differentiation variable
     */
    Transformation taylor(momentum) {
        use(Redberry) {
            { expr -> ((Differentiate[momentum.t] & momentum.t.eq(0.t)) >> expr) * momentum.t } as Transformation
        }
    }

    Map quarkoniaVertices = null

    /**
     * Calculating effective vertex g + g -> meson (scalar, axial, tensor)
     * @return map
     */
    Map effectiveQuarkoniaVertices() {
        if (quarkoniaVertices != null)
            return quarkoniaVertices;

        log 'Calculating effective vertex ...'
        use(Redberry) {
            def effVertices = [:]

            def factor = Factor[[FactorScalars: true, FactorizationEngine: wolframFactorTr]]

            // main formula (two diagrams)
            def A = 'cu[p1_m[fl]]*(V_aA*D[p1_m[fl] - k1_m, m[fl]]*V_bB + V_bB*D[p1_m[fl] - k2_m, m[fl]]*V_aA)*v[p2_m[fl]]'.t
            // basic simplifications
            A <<= FeynmanRules & spinSingletProjector['fl'] &
                    dTrace & ExpandAndEliminate & uTrace & momentums['fl'] &
                    ExpandAll[EliminateMetrics] & EliminateMetrics &
                    'p_m[fl]*p^m[fl] = (2*m[fl])**2'.t

            // Taylor expansion (scalar meson)
            def AScalar = ('q_i[fl] = q_i'.t.hold & taylor('q_i') & 'q_i = q_i[fl]'.t
                    & ExpandAll[EliminateMetrics] & EliminateMetrics
                    & 'p_m[fl]*p^m[fl] = (2*m[fl])**2'.t
                    & totalSpinProjector['scalar']
                    & 'p_i[fl] = k1_i + k2_i'.t
                    & ExpandAll[EliminateMetrics] & EliminateMetrics & Together
                    & factor
            ) >> A

            // Define
            effVertices['scalar'] = 'Ascalar_{aA bB}[fl, k1_m, k2_m]'.t.eq AScalar

            // Taylor expansion (axial meson)
            def AAxial = ('q_i[fl] = q_i'.t.hold & taylor('q_i') & 'q_i = q_i[fl]'.t
                    & ExpandAll[EliminateMetrics] & EliminateMetrics
                    & 'p_m[fl]*p^m[fl] = (2*m[fl])**2'.t
                    & totalSpinProjector['axial']
                    & Together
                    & ExpandAll[EliminateMetrics] & EliminateMetrics & Together
                    & factor
            ) >> A
            // Define
            effVertices['axial'] = 'Aaxial_{aA bB}[fl, k1_m, k2_m]'.t.eq AAxial

            // Taylor expansion (tensor meson)
            def ATensor = ('q_i[fl] = q_i'.t.hold & taylor('q_i') & 'q_i = q_i[fl]'.t
                    & ExpandAll[EliminateMetrics] & EliminateMetrics
                    & 'p_m[fl]*p^m[fl] = (2*m[fl])**2'.t
                    & totalSpinProjector['tensor']
                    & ExpandAll[EliminateMetrics] & EliminateMetrics
                    & 'p_m[fl]*p^m[fl] = (2*m[fl])**2'.t
                    & Together
                    & ExpandAll[EliminateMetrics] & EliminateMetrics & Together
                    & factor
            ) >> A
            // Define
            effVertices['tensor'] = 'Atensor_{aA bB}[fl, k1_m, k2_m]'.t.eq ATensor

            log '...done'
            return quarkoniaVertices = Collections.unmodifiableMap(effVertices)
        }
    }

    Transformation pairVertex = null

    /**
     * Effective vertex g + g -> q + qbar
     */
    Transformation effectivePairVertex() {
        if (pairVertex != null)
            return pairVertex

        use(Redberry) {
            //two diagrams
            def B = 'cu[p1_m[fl]]*(V_aA*D[p1_m[fl] - k1_m, m[fl]]*V_bB + V_bB*D[p1_m[fl] - k2_m, m[fl]]*V_aA)*v[p2_m[fl]]'.t

            // Simplifying
            B <<= FeynmanRules & ExpandAll[EliminateMetrics] & EliminateMetrics &
                    'p1_m[fl]*p1^m[fl] = m[fl]**2'.t & 'p2_m[fl]*p2^m[fl] = m[fl]**2'.t & Together

            return pairVertex = 'B_{aA bB}[fl, k1_m, k2_m]'.t.eq(B)
        }
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////// POLARIZATIONS ///////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////

    void setupPolarizationFactor() {
        use(Redberry) {
            if (projectCC)
                overallPolarizationFactor = 's*(16*mb**4 - 4*mb**2*s - 4*mb**2*t - 4*mb**2*u + t*u)'.t
            else
                overallPolarizationFactor = 's*(4*mc**4 - 4*mb**2*s - 4*mc**2*s + s**2 - 2*mc**2*t1 + s*t1 - 2*mc**2*t2 + s*t2 - 2*mc**2*u1 + s*u1 + t1*u1 + t2*u1 - 2*mc**2*u2 + s*u2 + t1*u2 + t2*u2)'.t
        }
    }

    Transformation setupPolarisations(def g1, def g2) {
        checkPol g1
        checkPol g2
        log 'Setting up polarizations'
        use(Redberry) {
            def cfs = Identity
            def eps1 = 'eps1_a = c1 * k1_a + c2 * k2_a + c3 * p_a[bottom]'.t
            def eps2 = 'eps2_a = c4 * e_abcd * k1^b * k2^c * p^d[bottom]'.t
            def epsPlus = 'eps_a[1] = eps1_a '.t << eps1 << eps2
            def epsMinus = 'eps_a[-1] = eps2_a'.t << eps1 << eps2

            if (projectCC) {
                cfs &= 'c1 = 4*mb**2 - t'.t
                cfs &= 'c2 = 4*mb**2 - u'.t
                cfs &= 'c3 = -s'.t
                cfs &= 'c4 = -2'.t
            } else {
                cfs &= 'c1 = -2*mc**2 + s + u1 + u2'.t
                cfs &= 'c2 = -2*mc**2 + s + t1 + t2'.t
                cfs &= 'c3 = -s'.t
                cfs &= 'c4 = -2'.t
            }
            epsPlus <<= cfs; epsMinus <<= cfs;

            def subs = Identity
            if (g1 != null)
                subs &= "eps1_a[h1] = eps_a[$g1]".t & epsPlus & epsMinus
            if (g2 != null)
                subs &= "eps2_a[h2] = eps_a[$g2]".t & epsPlus & epsMinus

            return subs
        }
    }

    def qqXYZCoeffs = [:]

    def calcQQXYZCoeffs(fl) {
        if (qqXYZCoeffs[fl] != null)
            return qqXYZCoeffs[fl]
        use(Redberry) {
            def var = fl[0] + 's'

            def eps1 = "eps1_a = ${var}1 * p_a[charm] + ${var}2 * p_a[bottom]".t
            def eps2 = "eps2_a = ${var}3 * p_a[charm] + ${var}4 * p_a[bottom] + ${var}5 * k1_a".t
            def eps0 = "eps0_a = ${var}6 * e_abcd * k1^b * p^c[charm] * p^d[bottom]".t

            if (!projectCC) {
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
            eq = (eps1 & eps2 & eps0 & ExpandAndEliminate & leviSimplify &
                    ExpandAndEliminate & mandelstam & massesSubs) >> eq
            def solverOptions = [ExternalSolver: [
                    Solver: 'Mathematica',
                    Path  : '/Applications/Mathematica.app/Contents/MacOS']]
            def solutions = Reduce(eq,
                    ["${var}1", "${var}2", "${var}3", "${var}4", "${var}5", "${var}6"].t,
                    solverOptions)
            assert solutions.size() != 0

            def cfs = solutions[0]
            cfs = (wolframFactorTr & wolframFactorSqrtTr) >> cfs
            return (qqXYZCoeffs[fl] = cfs)
        }
    }

    Map getXYZVector(def fl, def xyz) {
        use(Redberry) {
            def var = fl[0] + 's', epss = [:]

            def subs = Identity
            if (!projectCC)
                subs = 'p_a[charm] = p1_a[charm] + p2_a[charm]'.t.hold
            epss['x'] = "eps1_a = ${var}1 * p_a[charm] + ${var}2 * p_a[bottom]".t << subs
            epss['y'] = "eps2_a = ${var}3 * p_a[charm] + ${var}4 * p_a[bottom] + ${var}5 * k1_a".t << subs
            epss['z'] = "eps0_a = ${var}6 * e_abcd * k1^b * p^c[charm] * p^d[bottom]".t << subs


            def coeffs = calcQQXYZCoeffs(fl)
            def eps = epss[xyz][1] << coeffs
            eps <<= Together
            return ['den': Denominator >> eps, 'num': Numerator >> eps]
        }
    }

    Map setXYZ(def fl, def spin, def angular) {
        use(Redberry) {
            def vSpin = getXYZVector(fl, spin),
                vAngular = getXYZVector(fl, angular)
            return [den: vSpin['den'] * vAngular['den'],
                    tr : bind(["eps_i[hS[$fl]]": vSpin['num']]).hold &
                            bind(["eps_i[hL[$fl]]": vAngular['num']]).hold]
        }
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////// CALC ///////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////

    def spinorStructures, spinorStructuresVars, spinorSquares, conjugateSpinorL

    void setupSpinorStructures() {
        if (spinorStructures != null)
            return
        use(Redberry) {
            log 'Setting up multiplication table for spinor structures'
            def subs = []
            subs << 'cu[p1_m[charm]]*g_AB*v[p2_m[charm]]'.t
            subs << 'cu[p1_m[charm]]*T_A*T_B*v[p2_m[charm]]'.t
//            subs << 'cu[p1_m[charm]]*f_ABC*T^C*v[p2_m[charm]]'.t
//            subs << 'cu[p1_m[charm]]*d_ABC*T^C*v[p2_m[charm]]'.t

            subs << 'cu[p1_m[charm]]*g_AB*G5*v[p2_m[charm]]'.t
            subs << 'cu[p1_m[charm]]*T_A*T_B*G5*v[p2_m[charm]]'.t
//            subs << 'cu[p1_m[charm]]*f_ABC*T^C*G5*v[p2_m[charm]]'.t
//            subs << 'cu[p1_m[charm]]*d_ABC*T^C*G5*v[p2_m[charm]]'.t

            subs << 'cu[p1_m[charm]]*g_AB*G^i*v[p2_m[charm]]'.t
            subs << 'cu[p1_m[charm]]*T_A*T_B*G^i*v[p2_m[charm]]'.t
//            subs << 'cu[p1_m[charm]]*f_ABC*T^C*G^i*v[p2_m[charm]]'.t
//            subs << 'cu[p1_m[charm]]*d_ABC*T^C*G^i*v[p2_m[charm]]'.t

            subs << 'cu[p1_m[charm]]*g_AB*G^i*G5*v[p2_m[charm]]'.t
            subs << 'cu[p1_m[charm]]*T_A*T_B*G^i*G5*v[p2_m[charm]]'.t
//            subs << 'cu[p1_m[charm]]*f_ABC*T^C*G^i*G5*v[p2_m[charm]]'.t
//            subs << 'cu[p1_m[charm]]*d_ABC*T^C*G^i*G5*v[p2_m[charm]]'.t

            subs << 'cu[p1_m[charm]]*g_AB*G^i*G^j*v[p2_m[charm]]'.t
            subs << 'cu[p1_m[charm]]*T_A*T_B*G^i*G^j*v[p2_m[charm]]'.t
//            subs << 'cu[p1_m[charm]]*f_ABC*T^C*G^i*G^j*v[p2_m[charm]]'.t
//            subs << 'cu[p1_m[charm]]*d_ABC*T^C*G^i*G^j*v[p2_m[charm]]'.t

            subs << 'cu[p1_m[charm]]*g_AB*G^i*G^j*G5*v[p2_m[charm]]'.t
            subs << 'cu[p1_m[charm]]*T_A*T_B*G^i*G^j*G5*v[p2_m[charm]]'.t

//            subs << 'cu[p1_m[charm]]*g_AB*G^i*G^j*G^k*v[p2_m[charm]]'.t
//            subs << 'cu[p1_m[charm]]*T_A*T_B*G^i*G^j*G^k*v[p2_m[charm]]'.t
//
//            subs << 'cu[p1_m[charm]]*g_AB*G^i*G^j*G^k*G5*v[p2_m[charm]]'.t
//            subs << 'cu[p1_m[charm]]*T_A*T_B*G^i*G^j*G^k*G5*v[p2_m[charm]]'.t

            spinorStructuresVars = []
            conjugateSpinorL = Identity
            for (int i = 0; i < subs.size(); ++i) {
                def l = "L${i + 1}${subs[i].indices.free}".t
                subs[i] = subs[i].eq l
                spinorStructuresVars << l
                conjugateSpinorL &= "$l = c$l".t
            }

            spinorStructures = subs as Transformation
            spinorSquares = Identity

            def conjugate = Conjugate
            conjugate &= '{i -> a, j -> b, k -> c, A -> C, B -> D}'.mapping
            conjugate &= Reverse[Matrix1, Matrix2]
            conjugate &= conjugateSpinors
            conjugate &= 'G5 = -G5'.t

            for (int i = 0; i < subs.size(); ++i)
                for (int j = 0; j < subs.size(); ++j) {
                    def lhs = subs[i][1] * ((conjugate & conjugateSpinorL) >> subs[j][1])
                    def rhs = subs[i][0] * (conjugate >> subs[j][0])
                    rhs <<= epsSum & uTrace & mandelstam & dTraceSimplify &
                            fullSimplify & massesSubs & uSimplify & massesSubs & wFactor
                    spinorSquares &= lhs.eq(rhs)
                }
        }
    }

    Tensor calcAmplitude(Tensor amp, polarizations = Identity) {
        setupSpinorStructures()
        use(Redberry) {
            def qVertices = effectiveQuarkoniaVertices().values() as Transformation
            if (!projectCC)
                qVertices &= effectivePairVertex()

            amp <<= FeynmanRules & qVertices

            //processing denominator
            def den = Denominator >> amp
            den <<= ExpandAndEliminate & mandelstam & mandelstam & massesSubs & wFactor
            assert isSymbolic(den)

            //processing numerator
            def ls = LeviCivitaSimplify.minkowski[[OverallSimplifications: ExpandTensors[simplifyMetrics] & simplifyMetrics]]
            def fsE = simplifyMetrics & ExpandTensors[simplifyMetrics] & simplifyMetrics &
                    ls & ExpandTensors[simplifyMetrics] & simplifyMetrics & massesSubs

            def num = Numerator >> amp
            num <<= 'eps_i[hS[bottom]]*p^i[bottom] = 0'.t
            num <<= 'eps_i[hL[bottom]]*p^i[bottom] = 0'.t
            num <<= 'eps_i[h[bottom]]*p^i[bottom] = 0'.t
            num <<= 'eps_ij[h[bottom]]*p^i[bottom] = 0'.t
            num <<= polarizations & fsE & uTrace & EliminateMetrics & massesSubs

            println 'a'
            //reducing spinor structures
//            num <<= 'G_a*G_b*G_c = g_ab*G_c-g_ac*G_b+g_bc*G_a-I*e_abcd*G5*G^d'.t

            num <<= 'I*e_abcd*G^d = -G5*G_a*G_b*G_c + g_ab*G5*G_c-g_ac*G5*G_b+g_bc*G5*G_a'.t
//            num <<= 'T_A*T_B = 1/6*g_AB + 1/2*(I*f_ABC + d_ABC)*T^C'.t
            num <<= momentumConservation
            num <<= ExpandTensorsAndEliminate & simplifyMetrics & ls & spinorsSimplify & diracSimplify & ls &  massesSubs
            num <<= ExpandTensorsAndEliminate & simplifyMetrics & ls & spinorsSimplify & diracSimplify & ls &  massesSubs
            num <<= 'I*f_ABC*T^C = T_A*T_B - T_B*T_A'.t
            num <<= 'd_ABC*T^C = T_A*T_B + T_B*T_A - g_AB/3'.t

            println 'b'
            //instead of fullSimplifyE
            num <<= fsE


            println 'c'
            num <<= Expand >> 'eps_a[h[bottom]]*k1^a = -eps_a[h[bottom]]*(k2^a+p1^a[charm]+p2^a[charm])'.t
            num <<= Expand >> 'eps_ab[h[bottom]]*k1^a*k1^b = -eps_ab[h[bottom]]*(k2^a+p1^a[charm]+p2^a[charm])*(k2^b+p1^b[charm]+p2^b[charm])'.t
            num <<= Expand >> 'eps_ab[h[bottom]]*k1^a = -eps_ab[h[bottom]]*(k2^a+p1^a[charm]+p2^a[charm])'.t
            num <<= 'G_a*G_b*G_c = g_ab*G_c-g_ac*G_b+g_bc*G_a-I*e_abcd*G5*G^d'.t
            num <<= ExpandTensorsAndEliminate & simplifyMetrics & ls & spinorsSimplify & diracSimplify & ls &  massesSubs

            println 'd'
            num <<= fsE

            //replacing spinor structures
            num <<= spinorStructures
            num <<= Collect[*spinorStructuresVars, wFactor, [ExpandSymbolic: false]]
            log "Amplitude (numerator) info: ${info(num)}"
            //num.each { println it.dataSubProduct }
            return num / den
        }
    }

    Tensor calcProcess(Collection diagrams, polarizations = Identity) {
        use(Redberry) {
            log "Calculating amplitudes (${diagrams.size()})"
            SumBuilder sb = new SumBuilder()
            diagrams.eachWithIndex { e, i ->
                log "Calculating $i-th amplitude"
                sb << calcAmplitude(e, polarizations)
            }
            def M = sb.build()//wFactor >> sb.build()
            return overallPolarizationFactor * squareMatrixElement(M)
        }
    }

    Tensor squareMatrixElement(Tensor matrixElement) {
        squareMatrixElement(matrixElement, null)
    }

    public boolean doMapleFactorOnAmp2 = true

    Tensor squareMatrixElement(Tensor matrixElement, String spins) {
        use(Redberry) {
            if (matrixElement.class == Complex)
                return matrixElement**2
            spins = spins == null ? '' : "($spins)"
            log "Squaring matrix element $spins: ${info(matrixElement)}"

            def wrap = { expr -> expr.class == Sum ? expr as List : [expr] }
            def mElements = wrap(matrixElement)
            def conjugate = Conjugate & InvertIndices & conjugateSpinorL
            //def conjugate = Conjugate & InvertIndices
            //conjugate &= Reverse[Matrix1, Matrix2]
            //conjugate &= conjugateSpinors

            //def cMatrixElement = conjugate >> matrixElement
            //return ParallelExpand.parallelExpand(matrixElement, cMatrixElement, this.&calcProduct as Transformation, { l -> this.log(l) })

            def result = new SumBuilder()
            PTuples([mElements.size(), mElements.size()], 'squaring').each { i, j ->
                def part = mElements[i],
                    cPart = conjugate >> mElements[j]
                assert part.class == Product && cPart.class == Product

                def num_p = Numerator >> part,
                    cNum_p = Numerator >> cPart
                def num = wrap(num_p.class == Product ? num_p.dataSubProduct : num_p),
                    cNum = wrap(cNum_p.class == Product ? cNum_p.dataSubProduct : cNum_p)

                def numSb = new SumBuilder()
                PTuples([num.size(), cNum.size()], '  squaring term').each { i1, j1 ->
                    numSb << calcProduct(num[i1] * cNum[j1])
                }


                def overallNum = (num_p.class == Product ? num_p.indexlessSubProduct : '1'.t) * (cNum_p.class == Product ? cNum_p.indexlessSubProduct : 1.t) * numSb.build()
                log 'done term:'
                assert isSymbolic(overallNum)
                if (doMapleFactorOnAmp2)
                    overallNum <<= mapleFactorTr
                log(info(overallNum))

                def den = (Denominator >> part) * (Denominator >> cPart)
                result << (overallNum / den)
            }

            return result.build()
        }
    }

    Iterable<List<Integer>> PTuples(List bounds, pref = '') {
        def tuples = Tuples(bounds)
        return new Iterable<List<Integer>>() {
            @Override
            Iterator<List<Integer>> iterator() {
                return new IteratorWithProgress<List<Integer>>(tuples.iterator(),
                        bounds.inject(1, { a, b -> a * b }),
                        { p -> log "$pref $p %" })
            }
        }
    }

    Tensor calcProduct(Tensor t) {
        use(Redberry) {
            if (t.class != Product)
            //return (epsSum & uTrace & mandelstam & dTraceSimplify & fullSimplify & massesSubs & uSimplify) >> t
                return (epsSum & spinorSquares & fullSimplifyE & massesSubs & uSimplify & massesSubs) >> t
            Product part = t
            def indexless = part.indexlessSubProduct
            def tensor = part.dataSubProduct
            tensor <<= spinorSquares & epsSum & fullSimplifyE & massesSubs & uSimplify & massesSubs
            assert isSymbolic(tensor)
            def res = indexless * tensor
            return res
        }
    }

    void log(String msg, func = System.out.&println, printElapsed = true) {
        if (log) {
            if (!printElapsed)
                func msg
            else {
                int dur = (System.currentTimeMillis() - startTime) / 1000
                String elapsed = String.format("%d:%02d:%02d", (int) (dur / 3600), (int) (dur % 3600) / 60, (int) (dur % 60))
                func("[Redberry $elapsed] " + msg)
            }
        }
    }

    @Override
    void close() throws Exception {
        if (mathematicaKernel != null)
            mathematicaKernel.close()
        if (mapleEngine != null)
            mapleEngine.stop()
    }

    private static String os() {
        def os = System.getProperty('os.name').toLowerCase()
        if (os.contains('mac'))
            return 'mac'
        else if (os.contains('nux'))
            return 'linux'
        return null
    }

    static void checkPol(g) {
        if (g != null && g != 1 && g != -1)
            throw new IllegalArgumentException()
    }

    String tsts(Tensor t) {
        use(Redberry) {
            t <<= 's = 1'.t & 't1 = 2'.t & 't2 = 3'.t & 'u1 = 4'.t & 'u2 = 5'.t &
                    'g = 1'.t & 'mc = 1'.t & 'mb = 2'.t &
                    'm[charm] = 1'.t.hold & 'm[bottom] = 2'.t.hold
            return t
        }
    }
}
