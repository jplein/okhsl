package okhsl

import (
	"errors"
	"fmt"
	"math"
	"strconv"
)

var errOutOfBounds = errors.New("out of bounds")
var errInvalidHexString = errors.New("invalid hex string")

// Public data types

// HSL colorspace, normalized to values between 0 and 1

type HSLNormalized struct {
	H float64
	S float64
	L float64
}

func (h HSLNormalized) ToHSL() HSL {
	return HSL{
		H: h.H * 360,
		S: h.S * 100,
		L: h.L * 100,
	}
}

func (h HSLNormalized) Validate() error {
	if h.H < 0 || h.H > 1. {
		return getOutOfBoundsError("H", 0., 1., h.H)
	}

	if h.S < 0 || h.S > 1. {
		return getOutOfBoundsError("S", 0., 1., h.S)
	}

	if h.L < 0 || h.L > 1. {
		return getOutOfBoundsError("L", 0., 1., h.L)
	}

	return nil
}

// HSL colorspace which may feel more natural to humans:
// H - ranges from 0 to 360
// S - ranges from 0 to 100
// L - ranges from 0 to 100

type HSL struct {
	H float64
	S float64
	L float64
}

func (h HSL) ToHSLNormalized() HSLNormalized {
	return HSLNormalized{
		H: h.H / 360,
		S: h.S / 100,
		L: h.L / 100,
	}
}

func (h HSL) Validate() error {
	if h.H < 0 || h.H > 360. {
		return getOutOfBoundsError("H", 0., 360., h.H)
	}

	if h.S < 0 || h.S > 100. {
		return getOutOfBoundsError("S", 0., 100., h.S)
	}

	if h.L < 0 || h.L > 100. {
		return getOutOfBoundsError("L", 0., 100., h.L)
	}

	return nil
}

// RGB colorspace, normalized to values between 0 and 1

type RGBNormalized struct {
	R float64
	G float64
	B float64
}

func (r RGBNormalized) ToRGB() RGB {
	return RGB{
		R: r.R * 255,
		G: r.G * 255,
		B: r.B * 255,
	}
}

func (r RGBNormalized) Validate() error {
	if r.R < 0 || r.R > 1. {
		return getOutOfBoundsError("R", 0., 1., r.R)
	}

	if r.G < 0 || r.G > 1. {
		return getOutOfBoundsError("G", 0., 1., r.G)
	}

	if r.B < 0 || r.B > 1. {
		return getOutOfBoundsError("B", 0., 1., r.B)
	}

	return nil
}

// RGB colorspace which may feel more natural to humans:
// R - ranges from 0 to 255
// G - ranges from 0 to 255
// B - ranges from 0 to 255

type RGB struct {
	R float64
	G float64
	B float64
}

func (r RGB) ToRGBNormalized() RGBNormalized {
	return RGBNormalized{
		R: r.R / 255,
		G: r.G / 255,
		B: r.B / 255,
	}
}

func (r RGB) Validate() error {
	if r.R < 0 || r.R > 255. {
		return getOutOfBoundsError("R", 0., 255., r.R)
	}

	if r.G < 0 || r.G > 255. {
		return getOutOfBoundsError("G", 0., 255., r.G)
	}

	if r.B < 0 || r.B > 255. {
		return getOutOfBoundsError("B", 0., 255., r.B)
	}

	return nil
}

func (r RGB) ToHex() (string, error) {
	err := r.Validate()
	if err != nil {
		return "", err
	}

	rInt := int(math.Round(r.R))
	gInt := int(math.Round(r.G))
	bInt := int(math.Round(r.B))

	c := rInt<<16 + gInt<<8 + bInt

	return fmt.Sprintf("#%6x", c), nil
}

func (r *RGB) FromHex(hex string) error {
	if len(hex) != 7 {
		return fmt.Errorf("%w: %s: expected a string of 7 hex characters preceded by a # character", errInvalidHexString, hex)
	}

	if hex[0] != '#' {
		return fmt.Errorf("%w: %s: expected first character to be a #", errInvalidHexString, hex)
	}

	hexNum := hex[1:]

	c, err := strconv.ParseInt(hexNum, 16, 64)

	if err != nil {
		return fmt.Errorf("%w: %s: %w", errInvalidHexString, hex, err)
	}

	r.R = float64((c & 0xff0000) >> 16)
	r.G = float64((c & 0x00ff00) >> 8)
	r.B = float64((c & 0x0000ff) >> 0)

	return nil
}

// Public functions

func OKHSLToSRGBNormalized(hsl HSLNormalized) (RGBNormalized, error) {
	err := hsl.Validate()

	if err != nil {
		return RGBNormalized{}, err
	}

	h := hsl.H
	s := hsl.S
	l := hsl.L

	if l == 1.0 {
		return RGBNormalized{R: 1.0, G: 1.0, B: 1.0}, nil
	} else if l == 0.0 {
		return RGBNormalized{R: 0.0, G: 0.0, B: 0.0}, nil
	}

	a_ := math.Cos(2.0 * math.Pi * h)
	b_ := math.Sin(2.0 * math.Pi * h)
	L := toeInv(l)

	cs := getCs(L, a_, b_)
	C0 := cs.c0
	CMid := cs.cMid
	CMax := cs.cMax

	// Interpolate the three values for C so that:
	// At s=0: dC/ds = C_0, C=0
	// At s=0.8: C=C_mid
	// At s=1.0: C=C_max

	mid := 0.8
	midInv := 1.25

	var C, t, k0, k1, k2 float64

	if s < mid {
		t = midInv * s

		k1 = mid * C0
		k2 = (1.0 - k1/CMid)

		C = t * k1 / (1.0 - k2*t)
	} else {
		t = (s - mid) / (1 - mid)

		k0 = CMid
		k1 = (1.0 - mid) * CMid * CMid * midInv * midInv / C0
		k2 = (1.0 - (k1)/(CMax-CMid))

		C = k0 + t*k1/(1.0-k2*t)
	}

	rgb := oklabToLinearSRGB(lab{L, C * a_, C * b_})
	rgb_ := RGBNormalized{
		R: srgbTransferFunction(rgb.R),
		G: srgbTransferFunction(rgb.G),
		B: srgbTransferFunction(rgb.B),
	}

	return rgb_, nil
}

func SRGBToOKHSLNormalized(rgb RGBNormalized) (HSLNormalized, error) {
	err := rgb.Validate()
	if err != nil {
		return HSLNormalized{}, err
	}

	lab := linearSGBToOKLAB(RGBNormalized{
		R: srgbTransferFunctionInv(rgb.R),
		G: srgbTransferFunctionInv(rgb.G),
		B: srgbTransferFunctionInv(rgb.B),
	})

	C := math.Sqrt(lab.a*lab.a + lab.b*lab.b)
	a_ := lab.a / C
	b_ := lab.b / C

	L := lab.l
	h := 0.5 + 0.5*math.Atan2(-lab.b, -lab.a)/math.Pi

	cs := getCs(L, a_, b_)
	C0 := cs.c0
	CMid := cs.cMid
	CMax := cs.cMax

	// Inverse of the interpolation in okhsl_to_srgb:

	mid := 0.8
	midInv := 1.25

	var s float64

	if C < CMid {
		k1 := mid * C0
		k2 := (1.0 - k1/CMid)

		t := C / (k1 + k2*C)
		s = t * mid
	} else {
		k0 := CMid
		k1 := (1.0 - mid) * CMid * CMid * midInv * midInv / C0
		k2 := (1.0 - (k1)/(CMax-CMid))

		t := (C - k0) / (k1 + k2*(C-k0))
		s = mid + (1.0-mid)*t
	}

	l := toe(L)

	return HSLNormalized{H: h, S: s, L: l}, nil
}

// Private data types

type lab struct {
	l float64
	a float64
	b float64
}
type cs struct {
	c0   float64
	cMid float64
	cMax float64
}
type lc struct {
	l float64
	c float64
}
type st struct {
	s float64
	t float64
}

// Private functions

func getOutOfBoundsError(desc string, min, max, value float64) error {
	return fmt.Errorf("%w: %s: expected value in range %f to %f but found %f", errOutOfBounds, desc, min, max, value)
}

// Finds the maximum saturation possible for a given hue that fits in sRGB
// Saturation here is defined as S = C/L
// a and b must be normalized so a^2 + b^2 == 1.
func computeMaxSaturation(a, b float64) float64 {
	// Max saturation will be when one of r, g or b goes below zero.
	// Select different coefficients depending on which component goes below zero first
	var k0, k1, k2, k3, k4, wl, wm, ws float64

	switch {
	case -1.88170328*a-0.80936493*b > 1:
		// Red component
		k0 = +1.19086277
		k1 = +1.76576728
		k2 = +0.59662641
		k3 = +0.75515197
		k4 = +0.56771245
		wl = +4.0767416621
		wm = -3.3077115913
		ws = +0.2309699292
	case 1.81444104*a-1.19445276*b > 1:
		// Green component
		k0 = +0.73956515
		k1 = -0.45954404
		k2 = +0.08285427
		k3 = +0.12541070
		k4 = +0.14503204
		wl = -1.2684380046
		wm = +2.6097574011
		ws = -0.3413193965
	default:
		// Blue component
		k0 = +1.35733652
		k1 = -0.00915799
		k2 = -1.15130210
		k3 = -0.50559606
		k4 = +0.00692167
		wl = -0.0041960863
		wm = -0.7034186147
		ws = +1.7076147010
	}

	// Approximate max saturation using a polynomial:
	S := k0 + k1*a + k2*b + k3*a*a + k4*a*b

	// Do one step Halley's method to get closer
	// this gives an error less than 10e6, except for some blue hues where the dS/dh is close to infinite
	// this should be sufficient for most applications, otherwise do two/three steps

	kl := +0.3963377774*a + 0.2158037573*b
	km := -0.1055613458*a - 0.0638541728*b
	ks := -0.0894841775*a - 1.2914855480*b

	{
		l_ := 1. + S*kl
		m_ := 1. + S*km
		s_ := 1. + S*ks

		l := l_ * l_ * l_
		m := m_ * m_ * m_
		s := s_ * s_ * s_

		ldS := 3. * kl * l_ * l_
		mdS := 3. * km * m_ * m_
		sdS := 3. * ks * s_ * s_

		ldS2 := 6. * kl * kl * l_
		mdS2 := 6. * km * km * m_
		sdS2 := 6. * ks * ks * s_

		f := wl*l + wm*m + ws*s
		f1 := wl*ldS + wm*mdS + ws*sdS
		f2 := wl*ldS2 + wm*mdS2 + ws*sdS2

		S -= f * f1 / (f1*f1 - 0.5*f*f2)
	}

	return S
}

func oklabToLinearSRGB(c lab) RGBNormalized {
	l_ := c.l + 0.3963377774*c.a + 0.2158037573*c.b
	m_ := c.l - 0.1055613458*c.a - 0.0638541728*c.b
	s_ := c.l - 0.0894841775*c.a - 1.2914855480*c.b

	l := l_ * l_ * l_
	m := m_ * m_ * m_
	s := s_ * s_ * s_

	rgb := RGBNormalized{
		R: +4.0767416621*l - 3.3077115913*m + 0.2309699292*s,
		G: -1.2684380046*l + 2.6097574011*m - 0.3413193965*s,
		B: -0.0041960863*l - 0.7034186147*m + 1.7076147010*s,
	}

	return rgb
}

func findCusp(a, b float64) lc {
	// First, find the maximum saturation (saturation S = C/L)
	SCusp := computeMaxSaturation(a, b)

	// Convert to linear sRGB to find the first point where at least one of r,g or b >= 1:
	rgbAtMax := oklabToLinearSRGB(lab{1, SCusp * a, SCusp * b})
	LCusp := math.Cbrt(1.0 / math.Max(math.Max(rgbAtMax.R, rgbAtMax.G), rgbAtMax.B))
	CCusp := LCusp * SCusp

	return lc{LCusp, CCusp}
}

// Finds intersection of the line defined by
// L = L0 * (1 - t) + t * L1;
// C = t * C1;
// a and b must be normalized so a^2 + b^2 == 1.
func findGamutIntersection(a, b, l1, c1, l0 float64) float64 {
	// Find the cusp of the gamut triangle
	cusp := findCusp(a, b)

	// Find the intersection for upper and lower half seprately
	var t float64
	if ((l1-l0)*cusp.c - (cusp.l-l0)*c1) <= 0. {
		// Lower half
		t = cusp.c * l0 / (c1*cusp.l + cusp.c*(l0-l1))
	} else {
		// Upper half
		// First intersect with triangle
		t = cusp.c * (l0 - 1.) / (c1*(cusp.l-1.) + cusp.c*(l0-l1))

		// Then one step Halley's method
		{
			dL := l1 - l0
			dC := c1

			kl := +0.3963377774*a + 0.2158037573*b
			km := -0.1055613458*a - 0.0638541728*b
			ks := -0.0894841775*a - 1.2914855480*b

			ldt := dL + dC*kl
			mdt := dL + dC*km
			sdt := dL + dC*ks

			// If higher accuracy is required, 2 or 3 iterations of the following block can be used:
			{
				L := l0*(1.-t) + t*l1
				C := t * c1

				l_ := L + C*kl
				m_ := L + C*km
				s_ := L + C*ks

				l := l_ * l_ * l_
				m := m_ * m_ * m_
				s := s_ * s_ * s_

				ldt := 3 * ldt * l_ * l_
				mdt := 3 * mdt * m_ * m_
				sdt := 3 * sdt * s_ * s_

				ldt2 := 6 * ldt * ldt * l_
				mdt2 := 6 * mdt * mdt * m_
				sdt2 := 6 * sdt * sdt * s_

				r := 4.0767416621*l - 3.3077115913*m + 0.2309699292*s - 1
				r1 := 4.0767416621*ldt - 3.3077115913*mdt + 0.2309699292*sdt
				r2 := 4.0767416621*ldt2 - 3.3077115913*mdt2 + 0.2309699292*sdt2

				ur := r1 / (r1*r1 - 0.5*r*r2)
				tr := -r * ur

				g := -1.2684380046*l + 2.6097574011*m - 0.3413193965*s - 1
				g1 := -1.2684380046*ldt + 2.6097574011*mdt - 0.3413193965*sdt
				g2 := -1.2684380046*ldt2 + 2.6097574011*mdt2 - 0.3413193965*sdt2

				ug := g1 / (g1*g1 - 0.5*g*g2)
				tg := -g * ug

				b := -0.0041960863*l - 0.7034186147*m + 1.7076147010*s - 1
				b1 := -0.0041960863*ldt - 0.7034186147*mdt + 1.7076147010*sdt
				b2 := -0.0041960863*ldt2 - 0.7034186147*mdt2 + 1.7076147010*sdt2

				ub := b1 / (b1*b1 - 0.5*b*b2)
				tb := -b * ub

				// Original ternary expressions:
				// t_r = u_r >= 0.f ? t_r : FLT_MAX;
				// t_g = u_g >= 0.f ? t_g : FLT_MAX;
				// t_b = u_b >= 0.f ? t_b : FLT_MAX;

				if ur < 0 {
					tr = math.MaxFloat64
				}
				if ug < 0 {
					tg = math.MaxFloat64
				}
				if ub < 0 {
					tb = math.MaxFloat64
				}

				t += math.Min(tr, math.Min(tg, tb))
			}
		}
	}

	return t
}

func toST(cusp lc) st {
	L := cusp.l
	C := cusp.c

	return st{C / L, C / (1 - L)}
}

func getSTMid(a_, b_ float64) st {
	S := 0.11516993 + 1./(+7.44778970+4.15901240*b_+
		a_*(-2.19557347+1.75198401*b_+
			a_*(-2.13704948-10.02301043*b_+
				a_*(-4.24894561+5.38770819*b_+4.69891013*a_))))

	T := 0.11239642 + 1./(+1.61320320-0.68124379*b_+
		a_*(+0.40370612+0.90148123*b_+
			a_*(-0.27087943+0.61223990*b_+
				a_*(+0.00299215-0.45399568*b_-0.14661872*a_))))

	return st{S, T}
}

func getCs(l, a_, b_ float64) cs {
	cusp := findCusp(a_, b_)

	CMax := findGamutIntersection(a_, b_, l, 1, l /* cusp */)
	STMax := toST(cusp)

	// Scale factor to compensate for the curved part of gamut shape:
	k := CMax / math.Min((l*STMax.s), (1-l)*STMax.t)

	var CMid float64
	{
		STMid := getSTMid(a_, b_)

		// Use a soft minimum function, instead of a sharp triangle shape to get a smooth value for chroma.
		Ca := l * STMid.s
		Cb := (1.0 - l) * STMid.t
		CMid = 0.9 * k * math.Sqrt(math.Sqrt(1.0/(1.0/(Ca*Ca*Ca*Ca)+1.0/(Cb*Cb*Cb*Cb))))
	}

	var C0 float64
	{
		// for C_0, the shape is independent of hue, so ST are constant. Values picked to roughly be the average values of ST.
		Ca := l * 0.4
		Cb := (1.0 - l) * 0.8

		// Use a soft minimum function, instead of a sharp triangle shape to get a smooth value for chroma.
		C0 = math.Sqrt(1.0 / (1.0/(Ca*Ca) + 1.0/(Cb*Cb)))
	}

	return cs{c0: C0, cMid: CMid, cMax: CMax}
}

func toeInv(x float64) float64 {
	k1 := 0.206
	k2 := 0.03
	k3 := (1 + k1) / (1 + k2)

	return (x*x + k1*x) / (k3 * (x + k2))
}

func srgbTransferFunction(a float64) float64 {
	if .0031308 >= a {
		return 12.92 * a
	}

	return 1.055*math.Pow(a, .4166666666666667) - .055
}

func srgbTransferFunctionInv(a float64) float64 {
	if .04045 < a {
		return math.Pow((a+.055)/1.055, 2.4)
	} else {
		return a / 12.92
	}
}

func linearSGBToOKLAB(c RGBNormalized) lab {
	l := 0.4122214708*c.R + 0.5363325363*c.G + 0.0514459929*c.B
	m := 0.2119034982*c.R + 0.6806995451*c.G + 0.1073969566*c.B
	s := 0.0883024619*c.R + 0.2817188376*c.G + 0.6299787005*c.B

	l_ := math.Cbrt(l)
	m_ := math.Cbrt(m)
	s_ := math.Cbrt(s)

	return lab{
		l: 0.2104542553*l_ + 0.7936177850*m_ - 0.0040720468*s_,
		a: 1.9779984951*l_ - 2.4285922050*m_ + 0.4505937099*s_,
		b: 0.0259040371*l_ + 0.7827717662*m_ - 0.8086757660*s_,
	}
}

// toe function for L_r.
func toe(x float64) float64 {
	const k1 float64 = 0.206

	const k2 float64 = 0.03

	const k3 float64 = (1. + k1) / (1. + k2)

	return 0.5 * (k3*x - k1 + math.Sqrt((k3*x-k1)*(k3*x-k1)+4*k2*k3*x))
}
