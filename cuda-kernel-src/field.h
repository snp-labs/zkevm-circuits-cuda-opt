

typedef unsigned __int128 u128;

static_assert(sizeof(size_t) == 8);
static_assert(sizeof(ulong) == 8);
static_assert(sizeof(ulong2) == 16);
static_assert(sizeof(ulong4) == 32);
static_assert(sizeof(u128) == 16);

template <ulong modolus_x, ulong modolus_y, ulong modolus_z, ulong modolus_w, ulong inv>
class field
{
private:
    ulong4 self;

    /// Compute a + b + carry, returning the result and the new carry over.
    inline __device__ ulong2 adc(ulong a, ulong b, ulong carry) const
    {
        const u128 temp_a = a;
        const u128 temp_b = b;
        const u128 temp_carry = carry;

        const u128 result = temp_a + temp_b + temp_carry;

        return make_ulong2((ulong)result, (ulong)(result >> 64));
    }

    /// Compute a - (b + borrow), returning the result and the new borrow.
    inline __device__ ulong2 sbb(ulong a, ulong b, ulong borrow) const
    {
        const u128 temp_a = a;
        const u128 temp_b = b;
        const u128 temp_borrow = borrow;

        const u128 result = temp_a - (temp_b + (temp_borrow >> 63));

        return make_ulong2((ulong)result, (ulong)(result >> 64));
    }

    /// Compute a + (b * c) + carry, returning the result and the new carry over.
    inline __device__ ulong2 mac(ulong a, ulong b, ulong c, ulong carry) const
    {
        const u128 temp_a = a;
        const u128 temp_b = b;
        const u128 temp_c = c;
        const u128 temp_carry = carry;

        const u128 result = temp_a + (temp_b * temp_c) + temp_carry;

        return make_ulong2((ulong)result, (ulong)(result >> 64));
    }

    // The Montgomery reduction here is based on Algorithm 14.32 in
    // Handbook of Applied Cryptography
    // <http://cacr.uwaterloo.ca/hac/about/chap14.pdf>.
    __device__ ulong4 montgomery_reduce(ulong _r0, ulong _r1, ulong _r2, ulong _r3, ulong _r4, ulong _r5, ulong _r6, ulong _r7) const
    {
        ulong r1, r2, r3, r4, r5, r6, r7, d0, d1, d2, d3, carry2, borrow;

        const ulong4 MODULUS = make_ulong4(
            modolus_x,
            modolus_y,
            modolus_z,
            modolus_w);

        const ulong INV = inv;

        {
            const u128 temp_k = _r0;
            const ulong k = ((temp_k * INV) << 64) >> 64;
            const ulong2 temp0 = mac(_r0, k, MODULUS.x, 0);
            const ulong2 temp1 = mac(_r1, k, MODULUS.y, temp0.y);
            const ulong2 temp2 = mac(_r2, k, MODULUS.z, temp1.y);
            const ulong2 temp3 = mac(_r3, k, MODULUS.w, temp2.y);
            const ulong2 temp4 = adc(_r4, 0, temp3.y);

            r1 = temp1.x;
            r2 = temp2.x;
            r3 = temp3.x;
            r4 = temp4.x;
            carry2 = temp4.y;
        }

        {
            const u128 temp_k = r1;
            const ulong k = ((temp_k * INV) << 64) >> 64;
            const ulong2 temp0 = mac(r1, k, MODULUS.x, 0);
            const ulong2 temp1 = mac(r2, k, MODULUS.y, temp0.y);
            const ulong2 temp2 = mac(r3, k, MODULUS.z, temp1.y);
            const ulong2 temp3 = mac(r4, k, MODULUS.w, temp2.y);
            const ulong2 temp4 = adc(_r5, carry2, temp3.y);

            r2 = temp1.x;
            r3 = temp2.x;
            r4 = temp3.x;
            r5 = temp4.x;
            carry2 = temp4.y;
        }

        {
            const u128 temp_k = r2;
            const ulong k = ((temp_k * INV) << 64) >> 64;
            const ulong2 temp0 = mac(r2, k, MODULUS.x, 0);
            const ulong2 temp1 = mac(r3, k, MODULUS.y, temp0.y);
            const ulong2 temp2 = mac(r4, k, MODULUS.z, temp1.y);
            const ulong2 temp3 = mac(r5, k, MODULUS.w, temp2.y);
            const ulong2 temp4 = adc(_r6, carry2, temp3.y);

            r3 = temp1.x;
            r4 = temp2.x;
            r5 = temp3.x;
            r6 = temp4.x;
            carry2 = temp4.y;
        }

        {
            const u128 temp_k = r3;
            const ulong k = ((temp_k * INV) << 64) >> 64;
            const ulong2 temp0 = mac(r3, k, MODULUS.x, 0);
            const ulong2 temp1 = mac(r4, k, MODULUS.y, temp0.y);
            const ulong2 temp2 = mac(r5, k, MODULUS.z, temp1.y);
            const ulong2 temp3 = mac(r6, k, MODULUS.w, temp2.y);
            const ulong2 temp4 = adc(_r7, carry2, temp3.y);

            r4 = temp1.x;
            r5 = temp2.x;
            r6 = temp3.x;
            r7 = temp4.x;
            carry2 = temp4.y;
        }

        // Result may be within MODULUS of the correct value
        {
            const ulong2 temp0 = sbb(r4, MODULUS.x, 0);
            const ulong2 temp1 = sbb(r5, MODULUS.y, temp0.y);
            const ulong2 temp2 = sbb(r6, MODULUS.z, temp1.y);
            const ulong2 temp3 = sbb(r7, MODULUS.w, temp2.y);
            const ulong2 temp4 = sbb(carry2, 0, temp3.y);

            d0 = temp0.x;
            d1 = temp1.x;
            d2 = temp2.x;
            d3 = temp3.x;
            borrow = temp4.y;
        }

        {
            const ulong2 temp0 = adc(d0, MODULUS.x & borrow, 0);
            const ulong2 temp1 = adc(d1, MODULUS.y & borrow, temp0.y);
            const ulong2 temp2 = adc(d2, MODULUS.z & borrow, temp1.y);
            const ulong2 temp3 = adc(d3, MODULUS.w & borrow, temp2.y);

            d0 = temp0.x;
            d1 = temp1.x;
            d2 = temp2.x;
            d3 = temp3.x;
        }

        return make_ulong4(d0, d1, d2, d3);
    }

    /// Subtracts `rhs` from `self`, returning the result.
    __device__ ulong4 sub(const ulong4 self, const ulong4 rhs) const
    {

        const ulong4 MODULUS = make_ulong4(
            modolus_x,
            modolus_y,
            modolus_z,
            modolus_w);

        const ulong2 temp0 = sbb(self.x, rhs.x, 0);
        const ulong2 temp1 = sbb(self.y, rhs.y, temp0.y);
        const ulong2 temp2 = sbb(self.z, rhs.z, temp1.y);
        const ulong2 temp3 = sbb(self.w, rhs.w, temp2.y);

        // If underflow occurred on the final limb, borrow = 0xfff...fff, otherwise
        // borrow = 0x000...000. Thus, we use it as a mask to conditionally add the
        // modulus.
        const ulong2 temp4 = adc(temp0.x, MODULUS.x & temp3.y, 0);
        const ulong2 temp5 = adc(temp1.x, MODULUS.y & temp3.y, temp4.y);
        const ulong2 temp6 = adc(temp2.x, MODULUS.z & temp3.y, temp5.y);
        const ulong2 temp7 = adc(temp3.x, MODULUS.w & temp3.y, temp6.y);

        return make_ulong4(temp4.x, temp5.x, temp6.x, temp7.x);
    }

    /// Adds `rhs` to `self`, returning the result.
    __device__ ulong4 add(const ulong4 self, const ulong4 rhs) const
    {

        const ulong4 MODULUS = make_ulong4(
            modolus_x,
            modolus_y,
            modolus_z,
            modolus_w);

        const ulong2 temp0 = adc(self.x, rhs.x, 0);
        const ulong2 temp1 = adc(self.y, rhs.y, temp0.y);
        const ulong2 temp2 = adc(self.z, rhs.z, temp1.y);
        const ulong2 temp3 = adc(self.w, rhs.w, temp2.y);

        /// Attempt to subtract the modulus, to ensure the value
        /// is smaller than the modulus.
        return sub(make_ulong4(temp0.x, temp1.x, temp2.x, temp3.x), MODULUS);
    }

    //
    // Multiplies `rhs` by `self`, returning the result.
    // Schoolbook multiplication
    //
    __device__ ulong4 mul(const ulong4 self, const ulong4 rhs) const
    {

        ulong r0, r1, r2, r3, r4, r5, r6, r7;

        {
            const ulong2 temp0 = mac(0, self.x, rhs.x, 0);
            const ulong2 temp1 = mac(0, self.x, rhs.y, temp0.y);
            const ulong2 temp2 = mac(0, self.x, rhs.z, temp1.y);
            const ulong2 temp3 = mac(0, self.x, rhs.w, temp2.y);

            r0 = temp0.x;
            r1 = temp1.x;
            r2 = temp2.x;
            r3 = temp3.x;
            r4 = temp3.y;
        }

        {
            const ulong2 temp0 = mac(r1, self.y, rhs.x, 0);
            const ulong2 temp1 = mac(r2, self.y, rhs.y, temp0.y);
            const ulong2 temp2 = mac(r3, self.y, rhs.z, temp1.y);
            const ulong2 temp3 = mac(r4, self.y, rhs.w, temp2.y);

            r1 = temp0.x;
            r2 = temp1.x;
            r3 = temp2.x;
            r4 = temp3.x;
            r5 = temp3.y;
        }

        {
            const ulong2 temp0 = mac(r2, self.z, rhs.x, 0);
            const ulong2 temp1 = mac(r3, self.z, rhs.y, temp0.y);
            const ulong2 temp2 = mac(r4, self.z, rhs.z, temp1.y);
            const ulong2 temp3 = mac(r5, self.z, rhs.w, temp2.y);

            r2 = temp0.x;
            r3 = temp1.x;
            r4 = temp2.x;
            r5 = temp3.x;
            r6 = temp3.y;
        }

        {
            const ulong2 temp0 = mac(r3, self.w, rhs.x, 0);
            const ulong2 temp1 = mac(r4, self.w, rhs.y, temp0.y);
            const ulong2 temp2 = mac(r5, self.w, rhs.z, temp1.y);
            const ulong2 temp3 = mac(r6, self.w, rhs.w, temp2.y);

            r3 = temp0.x;
            r4 = temp1.x;
            r5 = temp2.x;
            r6 = temp3.x;
            r7 = temp3.y;
        }

        return montgomery_reduce(r0, r1, r2, r3, r4, r5, r6, r7);
    }

    /// Squares this element.
    __device__ ulong4 square(const ulong4 self) const
    {

        ulong r0, r1, r2, r3, r4, r5, r6, r7;

        {
            const ulong2 temp0 = mac(0, self.x, self.y, 0);
            const ulong2 temp1 = mac(0, self.x, self.z, temp0.y);
            const ulong2 temp2 = mac(0, self.x, self.w, temp1.y);

            r1 = temp0.x;
            r2 = temp1.x;
            r3 = temp2.x;
            r4 = temp2.y;
        }

        {
            const ulong2 temp0 = mac(r3, self.y, self.z, 0);
            const ulong2 temp1 = mac(r4, self.y, self.w, temp0.y);

            r3 = temp0.x;
            r4 = temp1.x;
            r5 = temp1.y;
        }

        {
            const ulong2 temp0 = mac(r5, self.z, self.w, 0);

            r5 = temp0.x;
            r6 = temp0.y;
        }

        r7 = r6 >> 63;
        r6 = (r6 << 1) | (r5 >> 63);
        r5 = (r5 << 1) | (r4 >> 63);
        r4 = (r4 << 1) | (r3 >> 63);
        r3 = (r3 << 1) | (r2 >> 63);
        r2 = (r2 << 1) | (r1 >> 63);
        r1 = r1 << 1;

        {
            const ulong2 temp0 = mac(0, self.x, self.x, 0);
            const ulong2 temp1 = adc(0, r1, temp0.y);
            const ulong2 temp2 = mac(r2, self.y, self.y, temp1.y);
            const ulong2 temp3 = adc(0, r3, temp2.y);
            const ulong2 temp4 = mac(r4, self.z, self.z, temp3.y);
            const ulong2 temp5 = adc(0, r5, temp4.y);
            const ulong2 temp6 = mac(r6, self.w, self.w, temp5.y);
            const ulong2 temp7 = adc(0, r7, temp6.y);

            r0 = temp0.x;
            r1 = temp1.x;
            r2 = temp2.x;
            r3 = temp3.x;
            r4 = temp4.x;
            r5 = temp5.x;
            r6 = temp6.x;
            r7 = temp7.x;
        }

        return montgomery_reduce(r0, r1, r2, r3, r4, r5, r6, r7);
    }

    /// Squares this element.
    __device__ ulong4 neg(const ulong4 self) const
    {
        const ulong4 MODULUS = make_ulong4(
            modolus_x,
            modolus_y,
            modolus_z,
            modolus_w);

        // Subtract `self` from `MODULUS` to negate. Ignore the final
        // borrow because it cannot underflow; self is guaranteed to
        // be in the field.
        const ulong2 temp0 = sbb(MODULUS.x, self.x, 0);
        const ulong2 temp1 = sbb(MODULUS.y, self.y, temp0.y);
        const ulong2 temp2 = sbb(MODULUS.z, self.z, temp1.y);
        const ulong2 temp3 = sbb(MODULUS.w, self.w, temp2.y);

        // // `tmp` could be `MODULUS` if `self` was zero. Create a mask that is
        // // zero if `self` was zero, and `u64::max_value()` if self was nonzero.
        return (self.x | self.y | self.z | self.w) ? make_ulong4(temp0.x, temp1.x, temp2.x, temp3.x) : make_ulong4(0, 0, 0, 0);
    }

public:
    inline __device__ field()
    {
        self = make_ulong4(0, 0, 0, 0);
    }

    inline __device__ field(const ulong4 _self)
    {
        self = make_ulong4(_self.x, _self.y, _self.z, _self.w);
    }

    inline __device__ field operator+(const field &other) const
    {
        return field(add(self, other.self));
    }

    inline __device__ field operator-(const field &other) const
    {
        return field(sub(self, other.self));
    }

    inline __device__ field operator*(const field &other) const
    {
        return field(mul(self, other.self));
    }

    inline __device__ field negate() const
    {
        return field(neg(self));
    }

    inline __device__ field square() const
    {
        return field(square(self));
    }
};

typedef field<(ulong)0x43e1f593f0000001,
              (ulong)0x2833e84879b97091,
              (ulong)0xb85045b68181585d,
              (ulong)0x30644e72e131a029,
              (ulong)0xc2e1f593efffffff>
    Fr;

typedef field<(ulong)0x3c208c16d87cfd47,
              (ulong)0x97816a916871ca8d,
              (ulong)0xb85045b68181585d,
              (ulong)0x30644e72e131a029,
              (ulong)0x87d20782e4866389>
    Fq;

static_assert(sizeof(ulong4) == sizeof(Fr));
static_assert(sizeof(ulong4) == sizeof(Fq));
