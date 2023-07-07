package okhsl

import (
	"math"
)

// functions below are ported from here:
// https://github.com/bottosson/bottosson.github.io/blob/master/misc/colorpicker/colorconversion.js

// Exported functions

func SRGBToOKHSL(r, g, b float64) []float64 {
	lab := linear_srgb_to_oklab(
		srgb_transfer_function_inv(r/255),
		srgb_transfer_function_inv(g/255),
		srgb_transfer_function_inv(b/255),
	)

	C := math.Sqrt(lab[1]*lab[1] + lab[2]*lab[2])
	a_ := lab[1] / C
	b_ := lab[2] / C

	L := lab[0]
	h := 0.5 + 0.5*math.Atan2(-lab[2], -lab[1])/math.Pi

	Cs := get_Cs(L, a_, b_)
	C_0 := Cs[0]
	C_mid := Cs[1]
	C_max := Cs[2]

	var s float64
	if C < C_mid {
		var k_0 float64 = 0
		k_1 := 0.8 * C_0
		k_2 := (1 - k_1/C_mid)

		t := (C - k_0) / (k_1 + k_2*(C-k_0))
		s = t * 0.8
	} else {
		k_0 := C_mid
		k_1 := 0.2 * C_mid * C_mid * 1.25 * 1.25 / C_0
		k_2 := (1 - (k_1)/(C_max-C_mid))

		t := (C - k_0) / (k_1 + k_2*(C-k_0))
		s = 0.8 + 0.2*t
	}

	l := toe(L)
	return []float64{h, s, l}
}

func OKHSLToSRGB(h, s, l float64) []float64 {
	if l == 1 {
		return []float64{255, 255, 255}
	} else if l == 0 {
		return []float64{0, 0, 0}
	}

	a_ := math.Cos(2 * math.Pi * h)
	b_ := math.Sin(2 * math.Pi * h)
	L := toe_inv(l)

	Cs := get_Cs(L, a_, b_)
	C_0 := Cs[0]
	C_mid := Cs[1]
	C_max := Cs[2]

	var C, t, k_0, k_1, k_2 float64
	if s < 0.8 {
		t = 1.25 * s
		k_0 = 0
		k_1 = 0.8 * C_0
		k_2 = (1 - k_1/C_mid)
	} else {
		t = 5 * (s - 0.8)
		k_0 = C_mid
		k_1 = 0.2 * C_mid * C_mid * 1.25 * 1.25 / C_0
		k_2 = (1 - (k_1)/(C_max-C_mid))
	}

	C = k_0 + t*k_1/(1-k_2*t)

	// If we would only use one of the Cs:
	//C = s*C_0;
	//C = s*1.25*C_mid;
	//C = s*C_max;

	rgb := oklab_to_linear_srgb(L, C*a_, C*b_)
	r := 255 * srgb_transfer_function(rgb[0])
	g := 255 * srgb_transfer_function(rgb[1])
	b := 255 * srgb_transfer_function(rgb[2])

	return []float64{r, g, b}
}

// Internal functions

func toe_inv(x float64) float64 {
	k_1 := 0.206
	k_2 := 0.03
	k_3 := (1 + k_1) / (1 + k_2)
	return (x*x + k_1*x) / (k_3 * (x + k_2))
}

// Finds the maximum saturation possible for a given hue that fits in sRGB
// Saturation here is defined as S = C/L
// a and b must be normalized so a^2 + b^2 == 1
func compute_max_saturation(a, b float64) float64 {
	// Max saturation will be when one of r, g or b goes below zero.

	// Select different coefficients depending on which component goes below zero first
	var k0, k1, k2, k3, k4, wl, wm, ws float64

	if -1.88170328*a-0.80936493*b > 1 {
		// Red component
		k0 = +1.19086277
		k1 = +1.76576728
		k2 = +0.59662641
		k3 = +0.75515197
		k4 = +0.56771245
		wl = +4.0767416621
		wm = -3.3077115913
		ws = +0.2309699292
	} else if 1.81444104*a-1.19445276*b > 1 {
		// Green component
		k0 = +0.73956515
		k1 = -0.45954404
		k2 = +0.08285427
		k3 = +0.12541070
		k4 = +0.14503204
		wl = -1.2684380046
		wm = +2.6097574011
		ws = -0.3413193965
	} else {
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

	k_l := +0.3963377774*a + 0.2158037573*b
	k_m := -0.1055613458*a - 0.0638541728*b
	k_s := -0.0894841775*a - 1.2914855480*b

	{
		l_ := 1 + S*k_l
		m_ := 1 + S*k_m
		s_ := 1 + S*k_s

		l := l_ * l_ * l_
		m := m_ * m_ * m_
		s := s_ * s_ * s_

		l_dS := 3 * k_l * l_ * l_
		m_dS := 3 * k_m * m_ * m_
		s_dS := 3 * k_s * s_ * s_

		l_dS2 := 6 * k_l * k_l * l_
		m_dS2 := 6 * k_m * k_m * m_
		s_dS2 := 6 * k_s * k_s * s_

		f := wl*l + wm*m + ws*s
		f1 := wl*l_dS + wm*m_dS + ws*s_dS
		f2 := wl*l_dS2 + wm*m_dS2 + ws*s_dS2

		S = S - f*f1/(f1*f1-0.5*f*f2)
	}

	return S
}

func oklab_to_linear_srgb(L, a, b float64) []float64 {
	l_ := L + 0.3963377774*a + 0.2158037573*b
	m_ := L - 0.1055613458*a - 0.0638541728*b
	s_ := L - 0.0894841775*a - 1.2914855480*b

	l := l_ * l_ * l_
	m := m_ * m_ * m_
	s := s_ * s_ * s_

	return []float64{
		(+4.0767416621*l - 3.3077115913*m + 0.2309699292*s),
		(-1.2684380046*l + 2.6097574011*m - 0.3413193965*s),
		(-0.0041960863*l - 0.7034186147*m + 1.7076147010*s),
	}
}

func find_cusp(a, b float64) []float64 {
	// First, find the maximum saturation (saturation S = C/L)
	S_cusp := compute_max_saturation(a, b)

	// Convert to linear sRGB to find the first point where at least one of r,g or b >= 1:
	rgb_at_max := oklab_to_linear_srgb(1, S_cusp*a, S_cusp*b)
	L_cusp := math.Cbrt(1 / math.Max(math.Max(rgb_at_max[0], rgb_at_max[1]), rgb_at_max[2]))
	C_cusp := L_cusp * S_cusp

	return []float64{L_cusp, C_cusp}
}

// Finds intersection of the line defined by
// L = L0 * (1 - t) + t * L1;
// C = t * C1;
// a and b must be normalized so a^2 + b^2 == 1
func find_gamut_intersection(a, b, L1, C1, L0 float64, cuspPtr *[]float64) float64 {
	var cusp []float64
	if cuspPtr == nil {
		// Find the cusp of the gamut triangle
		cusp = find_cusp(a, b)
	} else {
		cusp = *cuspPtr
	}

	// Find the intersection for upper and lower half seprately
	var t float64
	if ((L1-L0)*cusp[1] - (cusp[0]-L0)*C1) <= 0 {
		// Lower half

		t = cusp[1] * L0 / (C1*cusp[0] + cusp[1]*(L0-L1))
	} else {
		// Upper half

		// First intersect with triangle
		t = cusp[1] * (L0 - 1) / (C1*(cusp[0]-1) + cusp[1]*(L0-L1))

		// Then one step Halley's method
		{
			dL := L1 - L0
			dC := C1

			k_l := +0.3963377774*a + 0.2158037573*b
			k_m := -0.1055613458*a - 0.0638541728*b
			k_s := -0.0894841775*a - 1.2914855480*b

			l_dt := dL + dC*k_l
			m_dt := dL + dC*k_m
			s_dt := dL + dC*k_s

			// If higher accuracy is required, 2 or 3 iterations of the following block can be used:
			{
				L := L0*(1-t) + t*L1
				C := t * C1

				l_ := L + C*k_l
				m_ := L + C*k_m
				s_ := L + C*k_s

				l := l_ * l_ * l_
				m := m_ * m_ * m_
				s := s_ * s_ * s_

				ldt := 3 * l_dt * l_ * l_
				mdt := 3 * m_dt * m_ * m_
				sdt := 3 * s_dt * s_ * s_

				ldt2 := 6 * l_dt * l_dt * l_
				mdt2 := 6 * m_dt * m_dt * m_
				sdt2 := 6 * s_dt * s_dt * s_

				r := 4.0767416621*l - 3.3077115913*m + 0.2309699292*s - 1
				r1 := 4.0767416621*ldt - 3.3077115913*mdt + 0.2309699292*sdt
				r2 := 4.0767416621*ldt2 - 3.3077115913*mdt2 + 0.2309699292*sdt2

				u_r := r1 / (r1*r1 - 0.5*r*r2)
				t_r := -r * u_r

				g := -1.2684380046*l + 2.6097574011*m - 0.3413193965*s - 1
				g1 := -1.2684380046*ldt + 2.6097574011*mdt - 0.3413193965*sdt
				g2 := -1.2684380046*ldt2 + 2.6097574011*mdt2 - 0.3413193965*sdt2

				u_g := g1 / (g1*g1 - 0.5*g*g2)
				t_g := -g * u_g

				b := -0.0041960863*l - 0.7034186147*m + 1.7076147010*s - 1
				b1 := -0.0041960863*ldt - 0.7034186147*mdt + 1.7076147010*sdt
				b2 := -0.0041960863*ldt2 - 0.7034186147*mdt2 + 1.7076147010*sdt2

				u_b := b1 / (b1*b1 - 0.5*b*b2)
				t_b := -b * u_b

				if u_r < 0 {
					t_r = 10e5
				}

				if u_g < 0 {
					t_g = 10e5
				}

				if u_b < 0 {
					t_b = 10e5
				}

				t += math.Min(t_r, math.Min(t_g, t_b))
			}
		}
	}
	return t
}

func srgb_transfer_function(a float64) float64 {
	if .0031308 >= a {
		return 12.92 * a
	}
	return 1.055*math.Pow(a, .4166666666666667) - .055
}

func get_ST_max(a_, b_ float64, cuspPtr *[]float64) []float64 {
	var cusp []float64
	if cuspPtr == nil {
		cusp = find_cusp(a_, b_)
	} else {
		cusp = *cuspPtr
	}

	L := cusp[0]
	C := cusp[1]
	return []float64{C / L, C / (1 - L)}
}

func get_Cs(L, a_, b_ float64) []float64 {
	cusp := find_cusp(a_, b_)

	C_max := find_gamut_intersection(a_, b_, L, 1, L, &cusp)
	ST_max := get_ST_max(a_, b_, &cusp)

	S_mid := 0.11516993 + 1/(+7.44778970+4.15901240*b_+a_*(-2.19557347+1.75198401*b_+a_*(-2.13704948-10.02301043*b_+a_*(-4.24894561+5.38770819*b_+4.69891013*a_))))

	T_mid := 0.11239642 + 1/(+1.61320320-0.68124379*b_+a_*(+0.40370612+0.90148123*b_+a_*(-0.27087943+0.61223990*b_+a_*(+0.00299215-0.45399568*b_-0.14661872*a_))))
	k := C_max / math.Min((L*ST_max[0]), (1-L)*ST_max[1])

	var C_mid float64
	{
		C_a := L * S_mid
		C_b := (1 - L) * T_mid

		C_mid = 0.9 * k * math.Sqrt(math.Sqrt(1/(1/(C_a*C_a*C_a*C_a)+1/(C_b*C_b*C_b*C_b))))
	}

	var C_0 float64
	{
		C_a := L * 0.4
		C_b := (1 - L) * 0.8

		C_0 = math.Sqrt(1 / (1/(C_a*C_a) + 1/(C_b*C_b)))
	}

	return []float64{C_0, C_mid, C_max}
}

func linear_srgb_to_oklab(r, g, b float64) []float64 {
	l := 0.4122214708*r + 0.5363325363*g + 0.0514459929*b
	m := 0.2119034982*r + 0.6806995451*g + 0.1073969566*b
	s := 0.0883024619*r + 0.2817188376*g + 0.6299787005*b

	l_ := math.Cbrt(l)
	m_ := math.Cbrt(m)
	s_ := math.Cbrt(s)

	return []float64{
		0.2104542553*l_ + 0.7936177850*m_ - 0.0040720468*s_,
		1.9779984951*l_ - 2.4285922050*m_ + 0.4505937099*s_,
		0.0259040371*l_ + 0.7827717662*m_ - 0.8086757660*s_,
	}
}

func srgb_transfer_function_inv(a float64) float64 {
	if .04045 < a {
		return math.Pow((a+.055)/1.055, 2.4)
	} else {
		return a / 12.92
	}
}

func toe(x float64) float64 {
	k_1 := 0.206
	k_2 := 0.03
	k_3 := (1 + k_1) / (1 + k_2)

	return 0.5 * (k_3*x - k_1 + math.Sqrt((k_3*x-k_1)*(k_3*x-k_1)+4*k_2*k_3*x))
}

// don't need this one
func rgb_to_hsl(r float64, g float64, b float64) []float64 {
	r /= 255
	g /= 255
	b /= 255

	max := math.Max(r, math.Max(g, b))
	min := math.Min(r, math.Min(g, b))

	l := (max + min) / 2

	var h, s float64

	if max == min {
		h = 0
		s = 0
	} else {
		d := max - min
		if l > 0.5 {
			s = d / (2 - max - min)
		} else {
			s = d / (max + min)
		}

		switch max {
		case r:
			h = (g - b) / d
			if g < b {
				h += 6
			}
		case g:
			h = (b-r)/d + 2
		case b:
			h = (r-g)/d + 4
		}
		h /= 6
	}

	return []float64{h, s, l}
}
