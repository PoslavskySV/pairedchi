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
import static cc.redberry.core.tensor.Tensors.setAntiSymmetric
import static cc.redberry.core.tensor.Tensors.setSymmetric
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
               dTrace, dTraceSimplify, dSimplify, uTrace, uSimplify, leviSimplify, fullSimplify, fullSimplifyE,
               conjugateSpinors, momentumConservation

    /**
     * Polarisations
     */
    public Tensor overallPolarizationFactor

    //Mathematica related
    KernelLink mathematicaKernel
    public Transformation wFactor, wolframFactorTr, mSimplify, wolframSimplifyTr, wolframFactorSqrtTr

    //Maple related
    Engine mapleEngine;
    public Transformation mFactor, mapleFactorTr

    public Setup(boolean projectCC) {
        this(projectCC, false)
    }

    public Setup(boolean projectCC, boolean log) {
        this.startTime = System.currentTimeMillis()
        this.log = log
        this.projectCC = projectCC
        init()
    }

    void init() {//basic initialization
        use(Redberry) {
            log 'Setting up basic tensors'

            //Setting up matrix objects
            defineMatrices 'T_A', Matrix2.matrix,//unitary matrices
                    'G_a', 'G5', Matrix1.matrix, //gamma matrices
                    'v[p_a]', 'u[p_a]', Matrix1.vector, Matrix2.vector, //final quark wave function
                    'cu[p_a]', 'cv[p_a]', Matrix1.covector, Matrix2.covector, //final antiquark wave function
                    'V_iA', Matrix1.matrix, Matrix2.matrix, //quark-gluon vertex
                    'D[p_m, mass]', Matrix1.matrix //quark propagator

            //Levi-Civita, SU(N) symmetric and structure tensors
            Quiet { setAntiSymmetric 'e_abcd' }
            Quiet { setSymmetric 'd_ABC' }
            Quiet { setAntiSymmetric 'f_ABC' }
            //polarization tensor
            Quiet { setSymmetric 'eps_ab[h]' }

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
                totalSpinProjector['scalar'] &= "q_i[$fl]*epsS_j[$fl] = -J_ij[p_i[$fl], 2*m[$fl]]".t
                totalSpinProjector['axial'] &= "q_i[$fl]*epsS_j[$fl] = e_ijab * p^a[$fl] * eps^b[h[$fl]]".t
                totalSpinProjector['tensor'] &= "q_i[$fl]*epsS_j[$fl] = eps_ij[h[$fl]]".t

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
            else
                mandelstam = setMandelstam5([k1_m: '0', k2_m: '0', 'p1_m[charm]': 'm[charm]', 'p2_m[charm]': 'm[charm]', 'p_m[bottom]': '2*m[bottom]'])

            /** Simplifications of polarizations */
            def simplifyPolarizations = 'eps1^a[h1] * k1_a = 0'.t &
                    'eps2^a[h2] * k2_a = 0'.t &
                    'eps^a[h[charm]] * p_a[charm] = 0'.t &
                    'eps^a[h[bottom]] * p_a[bottom] = 0'.t &
                    'eps^ab[h[charm]] * p_a[charm] = 0'.t &
                    'eps^ab[h[bottom]] * p_a[bottom] = 0'.t

            //substituting masses for factor
            massesSubs = 'm[charm] = mc'.t.hold & 'm[bottom] = mb'.t.hold

            /** Full simplification */
            def simplifyMetrics = EliminateMetrics & simplifyPolarizations & mandelstam &
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

            dTraceSimplify = DiracTrace[[Gamma: 'G_a', Gamma5: 'G5', Simplifications: fullSimplify]]

            dSimplify = Identity
            dSimplify &= 'G_a*G^a = 4'.t & 'G_a*G_b*G^a = -2*G_b'.t & 'G_a*G_b*G_c*G^a = 4*g_bc'.t
            dSimplify &= 'G_a*G_b*eps^ab[h[charm]] = 0'.t & 'G_a*G_b*eps^ab[h[bottom]] = 0'.t
            def momentums
            if (projectCC)
                momentums = ['p_i[bottom]', 'p_i[charm]', 'k1_i', 'k2_i'].t
            else
                momentums = ['p_i[bottom]', 'p1_i[charm]', 'p2_i[charm]', 'k1_i', 'k2_i'].t
            momentums.each {
                def a = it
                def b = '{_i -> _j}'.mapping >> a
                def c = '{_i -> ^i}'.mapping >> a
                dSimplify &= mandelstam >> "G^i * G^j * $a * $b = $a * $c".t
            }

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

        log 'Calculating effective vertex'
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
            //TODO factor for others!!!

            // Define
            effVertices['scalar'] = 'Ascalar_{aA bB}[fl, k1_m, k2_m]'.t.eq AScalar

            // Taylor expansion (axial meson)
            def AAxial = ('q_i[fl] = q_i'.t.hold & taylor('q_i') & 'q_i = q_i[fl]'.t
                    & ExpandAll[EliminateMetrics] & EliminateMetrics
                    & 'p_m[fl]*p^m[fl] = (2*m[fl])**2'.t
                    & totalSpinProjector['axial']
                    & Together
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
            ) >> A
            // Define
            effVertices['tensor'] = 'Atensor_{aA bB}[fl, k1_m, k2_m]'.t.eq ATensor

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
//            def epsPlus = 'eps_a[1] = (eps1_a + I * eps2_a)/2**(1/2)'.t << eps1 << eps2
//            def epsMinus = 'eps_a[-1] = (eps1_a - I * eps2_a)/2**(1/2)'.t << eps1 << eps2

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

    ///////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////// CALC ///////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////

    def spinorStructures, spinorStructuresVars, spinorSquares

    void setupSpinorStructures() {
        if (spinorStructures != null)
            return
        use(Redberry) {
            log 'Setting up multiplication table for spinor structures'
            def subs = []
            subs << 'cu[p1_m[charm]]*g_AB*v[p2_m[charm]]'.t
            subs << 'cu[p1_m[charm]]*T_A*T_B*v[p2_m[charm]]'.t
            subs << 'cu[p1_m[charm]]*G^i*g_AB*v[p2_m[charm]]'.t
            subs << 'cu[p1_m[charm]]*G^i*G5*g_AB*v[p2_m[charm]]'.t
            subs << 'cu[p1_m[charm]]*T_A*T_B*G^i*v[p2_m[charm]]'.t
            subs << 'cu[p1_m[charm]]*T_A*T_B*G^i*G5*v[p2_m[charm]]'.t
            subs << 'cu[p1_m[charm]]*G^i*G^j*g_AB*v[p2_m[charm]]'.t
            subs << 'cu[p1_m[charm]]*G^i*G^j*G5*g_AB*v[p2_m[charm]]'.t
            subs << 'cu[p1_m[charm]]*T_A*T_B*G^i*G^j*v[p2_m[charm]]'.t
            subs << 'cu[p1_m[charm]]*T_A*T_B*G^i*G^j*G5*v[p2_m[charm]]'.t

            subs << 'cu[p1_m[charm]]*f_ABC*T^C*G^i*v[p2_m[charm]]'.t
            subs << 'cu[p1_m[charm]]*d_ABC*T^C*G^i*v[p2_m[charm]]'.t

            spinorStructuresVars = []
            for (int i = 0; i < subs.size(); ++i) {
                def l = "L${i + 1}${subs[i].indices.free}".t
                subs[i] = subs[i].eq l
                spinorStructuresVars << l
            }

            spinorStructures = subs as Transformation
            spinorSquares = Identity

            def conjugate = Conjugate
            conjugate &= '{i -> a, j -> b, A -> C, B -> D}'.mapping
            conjugate &= Reverse[Matrix1, Matrix2]
            conjugate &= conjugateSpinors

            for (int i = 0; i < subs.size(); ++i)
                for (int j = 0; j < subs.size(); ++j) {
                    def lhs = subs[i][1] * (conjugate >> subs[j][1])
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
            def num = Numerator >> amp
            num <<= polarizations & fullSimplify & uTrace & EliminateMetrics & massesSubs

            //reducing spinor structures
            num <<= dSimplify
            num <<= EliminateMetrics & massesSubs
            num <<= 'G_a*G_b*G_c = g_ab*G_c-g_ac*G_b+g_bc*G_a-I*e_abcd*G5*G^d'.t
            num <<= fullSimplifyE & EliminateMetrics & dSimplify & massesSubs
            num = Transformation.Util.applyUntilUnchanged(num, 'G5*G_a = -G_a*G5'.t)

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

    Tensor squareMatrixElement(Tensor matrixElement, String spins) {
        use(Redberry) {
            if (matrixElement.class == Complex)
                return matrixElement**2
            spins = spins == null ? '' : "($spins)"
            log "Squaring matrix element $spins: ${info(matrixElement)}"

            def wrap = { expr -> expr.class == Sum ? expr as List : [expr] }
            def mElements = wrap(matrixElement)
            def conjugate = Conjugate & InvertIndices
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
                //overallNum <<= mapleFactorTr
                log(info(overallNum))

                def den = (Denominator >> part) * (Denominator >> cPart)

//                    println 'num info'
//                    println info(num)
//                    StringBuilder sb = new StringBuilder()
//                    sb.append("r := ").append(num.toString(OutputFormat.Maple)).append(":")

//                    new File('/Users/poslavsky/Projects/redberry/redberry-pairedchi/output/res.maple') << sb.toString()
//                    new File('/Users/poslavsky/Projects/redberry/redberry-pairedchi/output/res.redberry') << num.toString(OutputFormat.Redberry)

//                    println info(num)

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

    def fileTensors = new File('/Users/poslavsky/Projects/redberry/redberry-pairedchi/output/tensors.redberry')
    def fileRedberry = new File('/Users/poslavsky/Projects/redberry/redberry-pairedchi/output/exprs.redberry')
    def fileMaple = new File('/Users/poslavsky/Projects/redberry/redberry-pairedchi/output/exprs.maple')
    def fileMathematica = new File('/Users/poslavsky/Projects/redberry/redberry-pairedchi/output/exprs.m')
    def counter = 0

    Tensor calcProduct(Tensor t) {
        use(Redberry) {
            if (t.class != Product)
            //return (epsSum & uTrace & mandelstam & dTraceSimplify & fullSimplify & massesSubs & uSimplify) >> t
                return (epsSum & spinorSquares & fullSimplifyE & massesSubs & uSimplify & massesSubs) >> t
            Product part = t
            def indexless = part.indexlessSubProduct
            def tensor = part.dataSubProduct

//            println '  '
//            println '  '
//            println '  '
//            tensor.each {
//                println it
//            }
//
//            fileTensors << tensor.toString(OutputFormat.Redberry)
//            fileTensors << '\n'

            //tensor <<= epsSum & uTrace & mandelstam & dTraceSimplify &
            //        fullSimplify & massesSubs & uSimplify & massesSubs
            tensor <<= spinorSquares & epsSum & fullSimplifyE & massesSubs & uSimplify & massesSubs

//            def var = "expr${counter++}"
//            StringBuilder sb = new StringBuilder()
//            sb.append(var).append(' = ').append(tensor).append(';\n')
//            fileRedberry << sb.toString()
//            sb = new StringBuilder()
//            sb.append(var).append(' = ').append(tensor.toString(WolframMathematica)).append(';\n')
//            fileMathematica << sb.toString()
//            sb = new StringBuilder()
//            sb.append(var).append(' := ').append(tensor.toString(Maple)).append(':\n')
//            fileMaple << sb.toString()

            assert isSymbolic(tensor)
            def res = indexless * tensor
//            println 'r done: w'
//            res <<= wolframFactorTr
//            println 'w done'
            return res
        }
    }

    void log(String msg) {
        if (log) {
            int dur = (System.currentTimeMillis() - startTime) / 1000
            String elapsed = String.format("%d:%02d:%02d", (int) (dur / 3600), (int) (dur % 3600) / 60, (int) (dur % 60))
            println("[Redberry $elapsed] " + msg)
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
}
