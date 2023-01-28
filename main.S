	.file	"main.cpp"
# GNU C++14 (Ubuntu 9.4.0-1ubuntu1~20.04.1) version 9.4.0 (x86_64-linux-gnu)
#	compiled by GNU C version 9.4.0, GMP version 6.2.0, MPFR version 4.0.2, MPC version 1.1.0, isl version isl-0.22.1-GMP

# GGC heuristics: --param ggc-min-expand=100 --param ggc-min-heapsize=131072
# options passed:  -imultiarch x86_64-linux-gnu -D_GNU_SOURCE main.cpp
# -mtune=generic -march=x86-64 -auxbase-strip main.S -fomit-frame-pointer
# -fverbose-asm -fasynchronous-unwind-tables -fstack-protector-strong
# -Wformat -Wformat-security -fstack-clash-protection -fcf-protection
# options enabled:  -fPIC -fPIE -faggressive-loop-optimizations
# -fassume-phsa -fasynchronous-unwind-tables -fauto-inc-dec -fcommon
# -fdelete-null-pointer-checks -fdwarf2-cfi-asm -fearly-inlining
# -feliminate-unused-debug-types -fexceptions -ffp-int-builtin-inexact
# -ffunction-cse -fgcse-lm -fgnu-runtime -fgnu-unique -fident
# -finline-atomics -fipa-stack-alignment -fira-hoist-pressure
# -fira-share-save-slots -fira-share-spill-slots -fivopts
# -fkeep-static-consts -fleading-underscore -flifetime-dse
# -flto-odr-type-merging -fmath-errno -fmerge-debug-strings
# -fomit-frame-pointer -fpeephole -fplt -fprefetch-loop-arrays
# -freg-struct-return -fsched-critical-path-heuristic
# -fsched-dep-count-heuristic -fsched-group-heuristic -fsched-interblock
# -fsched-last-insn-heuristic -fsched-rank-heuristic -fsched-spec
# -fsched-spec-insn-heuristic -fsched-stalled-insns-dep -fschedule-fusion
# -fsemantic-interposition -fshow-column -fshrink-wrap-separate
# -fsigned-zeros -fsplit-ivs-in-unroller -fssa-backprop
# -fstack-clash-protection -fstack-protector-strong -fstdarg-opt
# -fstrict-volatile-bitfields -fsync-libcalls -ftrapping-math -ftree-cselim
# -ftree-forwprop -ftree-loop-if-convert -ftree-loop-im -ftree-loop-ivcanon
# -ftree-loop-optimize -ftree-parallelize-loops= -ftree-phiprop
# -ftree-reassoc -ftree-scev-cprop -funit-at-a-time -funwind-tables
# -fverbose-asm -fzero-initialized-in-bss -m128bit-long-double -m64 -m80387
# -malign-stringops -mavx256-split-unaligned-load
# -mavx256-split-unaligned-store -mfancy-math-387 -mfp-ret-in-387 -mfxsr
# -mglibc -mieee-fp -mlong-double-80 -mmmx -mno-sse4 -mpush-args -mred-zone
# -msse -msse2 -mstv -mtls-direct-seg-refs -mvzeroupper

	.text
	.section	.rodata
	.type	_ZStL19piecewise_construct, @object
	.size	_ZStL19piecewise_construct, 1
_ZStL19piecewise_construct:
	.zero	1
	.section	.text._ZNKSt6chrono8durationIlSt5ratioILl1ELl1000000000EEE5countEv,"axG",@progbits,_ZNKSt6chrono8durationIlSt5ratioILl1ELl1000000000EEE5countEv,comdat
	.align 2
	.weak	_ZNKSt6chrono8durationIlSt5ratioILl1ELl1000000000EEE5countEv
	.type	_ZNKSt6chrono8durationIlSt5ratioILl1ELl1000000000EEE5countEv, @function
_ZNKSt6chrono8durationIlSt5ratioILl1ELl1000000000EEE5countEv:
.LFB1117:
	.cfi_startproc
	endbr64	
	movq	%rdi, -8(%rsp)	# this, this
# /usr/include/c++/9/chrono:347: 	{ return __r; }
	movq	-8(%rsp), %rax	# this, tmp84
	movq	(%rax), %rax	# this_2(D)->__r, _3
# /usr/include/c++/9/chrono:347: 	{ return __r; }
	ret	
	.cfi_endproc
.LFE1117:
	.size	_ZNKSt6chrono8durationIlSt5ratioILl1ELl1000000000EEE5countEv, .-_ZNKSt6chrono8durationIlSt5ratioILl1ELl1000000000EEE5countEv
	.section	.text._ZNSt6chrono8durationIlSt5ratioILl1ELl1000000000EEEC1IlvEERKT_,"axG",@progbits,_ZNSt6chrono8durationIlSt5ratioILl1ELl1000000000EEEC1IlvEERKT_,comdat
	.align 2
	.weak	_ZNSt6chrono8durationIlSt5ratioILl1ELl1000000000EEEC1IlvEERKT_
	.type	_ZNSt6chrono8durationIlSt5ratioILl1ELl1000000000EEEC1IlvEERKT_, @function
_ZNSt6chrono8durationIlSt5ratioILl1ELl1000000000EEEC1IlvEERKT_:
.LFB1122:
	.cfi_startproc
	endbr64	
	movq	%rdi, -8(%rsp)	# this, this
	movq	%rsi, -16(%rsp)	# __rep, __rep
# /usr/include/c++/9/chrono:332: 	  : __r(static_cast<rep>(__rep)) { }
	movq	-16(%rsp), %rax	# __rep, tmp83
	movq	(%rax), %rdx	# *__rep_5(D), _1
	movq	-8(%rsp), %rax	# this, tmp84
	movq	%rdx, (%rax)	# _1, this_3(D)->__r
# /usr/include/c++/9/chrono:332: 	  : __r(static_cast<rep>(__rep)) { }
	nop	
	ret	
	.cfi_endproc
.LFE1122:
	.size	_ZNSt6chrono8durationIlSt5ratioILl1ELl1000000000EEEC1IlvEERKT_, .-_ZNSt6chrono8durationIlSt5ratioILl1ELl1000000000EEEC1IlvEERKT_
	.section	.text._ZSt4sqrtf,"axG",@progbits,_ZSt4sqrtf,comdat
	.weak	_ZSt4sqrtf
	.type	_ZSt4sqrtf, @function
_ZSt4sqrtf:
.LFB1210:
	.cfi_startproc
	endbr64	
	subq	$24, %rsp	#,
	.cfi_def_cfa_offset 32
	movss	%xmm0, 12(%rsp)	# __x, __x
# /usr/include/c++/9/cmath:464:   { return __builtin_sqrtf(__x); }
	movl	12(%rsp), %eax	# __x, tmp84
	movd	%eax, %xmm0	# tmp84,
	call	sqrtf@PLT	#
# /usr/include/c++/9/cmath:464:   { return __builtin_sqrtf(__x); }
	addq	$24, %rsp	#,
	.cfi_def_cfa_offset 8
	ret	
	.cfi_endproc
.LFE1210:
	.size	_ZSt4sqrtf, .-_ZSt4sqrtf
	.section	.rodata
	.align 8
	.type	_ZL1N, @object
	.size	_ZL1N, 8
_ZL1N:
	.quad	48
	.text
	.globl	_Z5frandv
	.type	_Z5frandv, @function
_Z5frandv:
.LFB1369:
	.cfi_startproc
	endbr64	
	subq	$8, %rsp	#,
	.cfi_def_cfa_offset 16
# main.cpp:10:     return (float)std::rand() / RAND_MAX * 2 - 1;
	call	rand@PLT	#
# main.cpp:10:     return (float)std::rand() / RAND_MAX * 2 - 1;
	cvtsi2ssl	%eax, %xmm0	# _1, _2
# main.cpp:10:     return (float)std::rand() / RAND_MAX * 2 - 1;
	movss	.LC0(%rip), %xmm1	#, tmp88
	divss	%xmm1, %xmm0	# tmp88, _3
# main.cpp:10:     return (float)std::rand() / RAND_MAX * 2 - 1;
	addss	%xmm0, %xmm0	# _3, _4
# main.cpp:10:     return (float)std::rand() / RAND_MAX * 2 - 1;
	movss	.LC1(%rip), %xmm1	#, tmp89
	subss	%xmm1, %xmm0	# tmp89, _7
# main.cpp:11: }
	addq	$8, %rsp	#,
	.cfi_def_cfa_offset 8
	ret	
	.cfi_endproc
.LFE1369:
	.size	_Z5frandv, .-_Z5frandv
	.globl	stars
	.bss
	.align 32
	.type	stars, @object
	.size	stars, 1344
stars:
	.zero	1344
	.text
	.globl	_Z4initv
	.type	_Z4initv, @function
_Z4initv:
.LFB1370:
	.cfi_startproc
	endbr64	
	subq	$24, %rsp	#,
	.cfi_def_cfa_offset 32
# main.cpp:23:     for (size_t i = 0; i < N; i++) {
	movq	$0, 8(%rsp)	#, i
.L10:
# main.cpp:23:     for (size_t i = 0; i < N; i++) {
	cmpq	$47, 8(%rsp)	#, i
	setbe	%al	#, _1
	testb	%al, %al	# _1
	je	.L11	#,
# main.cpp:24:         stars.px[i] = frand();
	call	_Z5frandv	#
	movd	%xmm0, %eax	#, _2
# main.cpp:24:         stars.px[i] = frand();
	movq	8(%rsp), %rdx	# i, tmp91
	leaq	0(,%rdx,4), %rcx	#, tmp92
	leaq	stars(%rip), %rdx	#, tmp93
	movl	%eax, (%rcx,%rdx)	# _2, stars.px
# main.cpp:25:         stars.py[i] = frand();
	call	_Z5frandv	#
	movd	%xmm0, %eax	#, _3
# main.cpp:25:         stars.py[i] = frand();
	movq	8(%rsp), %rdx	# i, tmp95
	addq	$48, %rdx	#, tmp94
	leaq	0(,%rdx,4), %rcx	#, tmp96
	leaq	stars(%rip), %rdx	#, tmp97
	movl	%eax, (%rcx,%rdx)	# _3, stars.py
# main.cpp:26:         stars.pz[i] = frand();
	call	_Z5frandv	#
	movd	%xmm0, %eax	#, _4
# main.cpp:26:         stars.pz[i] = frand();
	movq	8(%rsp), %rdx	# i, tmp99
	addq	$96, %rdx	#, tmp98
	leaq	0(,%rdx,4), %rcx	#, tmp100
	leaq	stars(%rip), %rdx	#, tmp101
	movl	%eax, (%rcx,%rdx)	# _4, stars.pz
# main.cpp:27:         stars.vx[i] = frand();
	call	_Z5frandv	#
	movd	%xmm0, %eax	#, _5
# main.cpp:27:         stars.vx[i] = frand();
	movq	8(%rsp), %rdx	# i, tmp103
	addq	$144, %rdx	#, tmp102
	leaq	0(,%rdx,4), %rcx	#, tmp104
	leaq	stars(%rip), %rdx	#, tmp105
	movl	%eax, (%rcx,%rdx)	# _5, stars.vx
# main.cpp:28:         stars.vy[i] = frand();
	call	_Z5frandv	#
	movd	%xmm0, %eax	#, _6
# main.cpp:28:         stars.vy[i] = frand();
	movq	8(%rsp), %rdx	# i, tmp107
	addq	$192, %rdx	#, tmp106
	leaq	0(,%rdx,4), %rcx	#, tmp108
	leaq	stars(%rip), %rdx	#, tmp109
	movl	%eax, (%rcx,%rdx)	# _6, stars.vy
# main.cpp:29:         stars.vz[i] = frand();
	call	_Z5frandv	#
	movd	%xmm0, %eax	#, _7
# main.cpp:29:         stars.vz[i] = frand();
	movq	8(%rsp), %rdx	# i, tmp111
	addq	$240, %rdx	#, tmp110
	leaq	0(,%rdx,4), %rcx	#, tmp112
	leaq	stars(%rip), %rdx	#, tmp113
	movl	%eax, (%rcx,%rdx)	# _7, stars.vz
# main.cpp:30:         stars.mass[i] = frand() + 1.f;
	call	_Z5frandv	#
# main.cpp:30:         stars.mass[i] = frand() + 1.f;
	movss	.LC1(%rip), %xmm1	#, tmp114
	addss	%xmm1, %xmm0	# tmp114, _9
# main.cpp:30:         stars.mass[i] = frand() + 1.f;
	movq	8(%rsp), %rax	# i, tmp116
	addq	$288, %rax	#, tmp115
	leaq	0(,%rax,4), %rdx	#, tmp117
	leaq	stars(%rip), %rax	#, tmp118
	movss	%xmm0, (%rdx,%rax)	# _9, stars.mass
# main.cpp:23:     for (size_t i = 0; i < N; i++) {
	addq	$1, 8(%rsp)	#, i
	jmp	.L10	#
.L11:
# main.cpp:32: }
	nop	
	addq	$24, %rsp	#,
	.cfi_def_cfa_offset 8
	ret	
	.cfi_endproc
.LFE1370:
	.size	_Z4initv, .-_Z4initv
	.section	.rodata
	.align 4
	.type	_ZL1G, @object
	.size	_ZL1G, 4
_ZL1G:
	.long	981668463
	.align 4
	.type	_ZL3eps, @object
	.size	_ZL3eps, 4
_ZL3eps:
	.long	981668463
	.align 4
	.type	_ZL2dt, @object
	.size	_ZL2dt, 4
_ZL2dt:
	.long	1008981770
	.align 4
	.type	_ZL4eps2, @object
	.size	_ZL4eps2, 4
_ZL4eps2:
	.long	897988542
	.align 4
	.type	_ZL3dt2, @object
	.size	_ZL3dt2, 4
_ZL3dt2:
	.long	953267991
	.align 4
	.type	_ZL9G_plus_dt, @object
	.size	_ZL9G_plus_dt, 4
_ZL9G_plus_dt:
	.long	925353389
	.text
	.globl	_Z4stepv
	.type	_Z4stepv, @function
_Z4stepv:
.LFB1371:
	.cfi_startproc
	endbr64	
	subq	$56, %rsp	#,
	.cfi_def_cfa_offset 64
# main.cpp:42:     for (size_t i = 0; i < N; i++) {
	movq	$0, 24(%rsp)	#, i
.L16:
# main.cpp:42:     for (size_t i = 0; i < N; i++) {
	cmpq	$47, 24(%rsp)	#, i
	ja	.L13	#,
# main.cpp:44:         for (size_t j = 0; j < N; j++) {
	movq	$0, 32(%rsp)	#, j
.L15:
# main.cpp:44:         for (size_t j = 0; j < N; j++) {
	cmpq	$47, 32(%rsp)	#, j
	setbe	%al	#, _1
	testb	%al, %al	# _1
	je	.L14	#,
# main.cpp:45:             float dx = stars.px[j] - stars.px[i];
	movq	32(%rsp), %rax	# j, tmp126
	leaq	0(,%rax,4), %rdx	#, tmp127
	leaq	stars(%rip), %rax	#, tmp128
	movss	(%rdx,%rax), %xmm0	# stars.px, _2
# main.cpp:45:             float dx = stars.px[j] - stars.px[i];
	movq	24(%rsp), %rax	# i, tmp129
	leaq	0(,%rax,4), %rdx	#, tmp130
	leaq	stars(%rip), %rax	#, tmp131
	movss	(%rdx,%rax), %xmm1	# stars.px, _3
# main.cpp:45:             float dx = stars.px[j] - stars.px[i];
	subss	%xmm1, %xmm0	# _3, tmp132
	movss	%xmm0, 4(%rsp)	# tmp132, dx
# main.cpp:46:             float dy = stars.py[j] - stars.py[i];
	movq	32(%rsp), %rax	# j, tmp134
	addq	$48, %rax	#, tmp133
	leaq	0(,%rax,4), %rdx	#, tmp135
	leaq	stars(%rip), %rax	#, tmp136
	movss	(%rdx,%rax), %xmm0	# stars.py, _4
# main.cpp:46:             float dy = stars.py[j] - stars.py[i];
	movq	24(%rsp), %rax	# i, tmp138
	addq	$48, %rax	#, tmp137
	leaq	0(,%rax,4), %rdx	#, tmp139
	leaq	stars(%rip), %rax	#, tmp140
	movss	(%rdx,%rax), %xmm1	# stars.py, _5
# main.cpp:46:             float dy = stars.py[j] - stars.py[i];
	subss	%xmm1, %xmm0	# _5, tmp141
	movss	%xmm0, 8(%rsp)	# tmp141, dy
# main.cpp:47:             float dz = stars.pz[j] - stars.pz[i];
	movq	32(%rsp), %rax	# j, tmp143
	addq	$96, %rax	#, tmp142
	leaq	0(,%rax,4), %rdx	#, tmp144
	leaq	stars(%rip), %rax	#, tmp145
	movss	(%rdx,%rax), %xmm0	# stars.pz, _6
# main.cpp:47:             float dz = stars.pz[j] - stars.pz[i];
	movq	24(%rsp), %rax	# i, tmp147
	addq	$96, %rax	#, tmp146
	leaq	0(,%rax,4), %rdx	#, tmp148
	leaq	stars(%rip), %rax	#, tmp149
	movss	(%rdx,%rax), %xmm1	# stars.pz, _7
# main.cpp:47:             float dz = stars.pz[j] - stars.pz[i];
	subss	%xmm1, %xmm0	# _7, tmp150
	movss	%xmm0, 12(%rsp)	# tmp150, dz
# main.cpp:48:             float d2 = dx * dx + dy * dy + dz * dz + eps2;
	movss	4(%rsp), %xmm0	# dx, tmp151
	movaps	%xmm0, %xmm1	# tmp151, tmp151
	mulss	%xmm0, %xmm1	# tmp151, tmp151
# main.cpp:48:             float d2 = dx * dx + dy * dy + dz * dz + eps2;
	movss	8(%rsp), %xmm0	# dy, tmp152
	mulss	%xmm0, %xmm0	# tmp152, _9
# main.cpp:48:             float d2 = dx * dx + dy * dy + dz * dz + eps2;
	addss	%xmm0, %xmm1	# _9, _10
# main.cpp:48:             float d2 = dx * dx + dy * dy + dz * dz + eps2;
	movss	12(%rsp), %xmm0	# dz, tmp153
	mulss	%xmm0, %xmm0	# tmp153, _11
# main.cpp:48:             float d2 = dx * dx + dy * dy + dz * dz + eps2;
	addss	%xmm0, %xmm1	# _11, _12
# main.cpp:48:             float d2 = dx * dx + dy * dy + dz * dz + eps2;
	movss	.LC2(%rip), %xmm0	#, tmp155
	addss	%xmm1, %xmm0	# _12, tmp154
	movss	%xmm0, 16(%rsp)	# tmp154, d2
# main.cpp:49:             d2 *= std::sqrt(d2);
	movl	16(%rsp), %eax	# d2, tmp156
	movd	%eax, %xmm0	# tmp156,
	call	_ZSt4sqrtf	#
# main.cpp:49:             d2 *= std::sqrt(d2);
	movss	16(%rsp), %xmm1	# d2, tmp158
	mulss	%xmm1, %xmm0	# tmp158, tmp157
	movss	%xmm0, 16(%rsp)	# tmp157, d2
# main.cpp:50:             float inv_d2 = 1.f / d2;
	movss	.LC1(%rip), %xmm0	#, tmp160
	divss	16(%rsp), %xmm0	# d2, tmp159
	movss	%xmm0, 20(%rsp)	# tmp159, inv_d2
# main.cpp:51:             stars.vx[i] += dx * stars.mass[j] * G_plus_dt * inv_d2;
	movq	24(%rsp), %rax	# i, tmp162
	addq	$144, %rax	#, tmp161
	leaq	0(,%rax,4), %rdx	#, tmp163
	leaq	stars(%rip), %rax	#, tmp164
	movss	(%rdx,%rax), %xmm1	# stars.vx, _13
# main.cpp:51:             stars.vx[i] += dx * stars.mass[j] * G_plus_dt * inv_d2;
	movq	32(%rsp), %rax	# j, tmp166
	addq	$288, %rax	#, tmp165
	leaq	0(,%rax,4), %rdx	#, tmp167
	leaq	stars(%rip), %rax	#, tmp168
	movss	(%rdx,%rax), %xmm0	# stars.mass, _14
# main.cpp:51:             stars.vx[i] += dx * stars.mass[j] * G_plus_dt * inv_d2;
	movaps	%xmm0, %xmm2	# _14, _14
	mulss	4(%rsp), %xmm2	# dx, _14
# main.cpp:51:             stars.vx[i] += dx * stars.mass[j] * G_plus_dt * inv_d2;
	movss	.LC3(%rip), %xmm0	#, tmp169
	mulss	%xmm2, %xmm0	# _15, _16
# main.cpp:51:             stars.vx[i] += dx * stars.mass[j] * G_plus_dt * inv_d2;
	mulss	20(%rsp), %xmm0	# inv_d2, _17
# main.cpp:51:             stars.vx[i] += dx * stars.mass[j] * G_plus_dt * inv_d2;
	addss	%xmm1, %xmm0	# _13, _18
	movq	24(%rsp), %rax	# i, tmp171
	addq	$144, %rax	#, tmp170
	leaq	0(,%rax,4), %rdx	#, tmp172
	leaq	stars(%rip), %rax	#, tmp173
	movss	%xmm0, (%rdx,%rax)	# _18, stars.vx
# main.cpp:52:             stars.vy[i] += dy * stars.mass[j] * G_plus_dt * inv_d2;
	movq	24(%rsp), %rax	# i, tmp175
	addq	$192, %rax	#, tmp174
	leaq	0(,%rax,4), %rdx	#, tmp176
	leaq	stars(%rip), %rax	#, tmp177
	movss	(%rdx,%rax), %xmm1	# stars.vy, _19
# main.cpp:52:             stars.vy[i] += dy * stars.mass[j] * G_plus_dt * inv_d2;
	movq	32(%rsp), %rax	# j, tmp179
	addq	$288, %rax	#, tmp178
	leaq	0(,%rax,4), %rdx	#, tmp180
	leaq	stars(%rip), %rax	#, tmp181
	movss	(%rdx,%rax), %xmm0	# stars.mass, _20
# main.cpp:52:             stars.vy[i] += dy * stars.mass[j] * G_plus_dt * inv_d2;
	movaps	%xmm0, %xmm2	# _20, _20
	mulss	8(%rsp), %xmm2	# dy, _20
# main.cpp:52:             stars.vy[i] += dy * stars.mass[j] * G_plus_dt * inv_d2;
	movss	.LC3(%rip), %xmm0	#, tmp182
	mulss	%xmm2, %xmm0	# _21, _22
# main.cpp:52:             stars.vy[i] += dy * stars.mass[j] * G_plus_dt * inv_d2;
	mulss	20(%rsp), %xmm0	# inv_d2, _23
# main.cpp:52:             stars.vy[i] += dy * stars.mass[j] * G_plus_dt * inv_d2;
	addss	%xmm1, %xmm0	# _19, _24
	movq	24(%rsp), %rax	# i, tmp184
	addq	$192, %rax	#, tmp183
	leaq	0(,%rax,4), %rdx	#, tmp185
	leaq	stars(%rip), %rax	#, tmp186
	movss	%xmm0, (%rdx,%rax)	# _24, stars.vy
# main.cpp:53:             stars.vz[i] += dz * stars.mass[j] * G_plus_dt * inv_d2;
	movq	24(%rsp), %rax	# i, tmp188
	addq	$240, %rax	#, tmp187
	leaq	0(,%rax,4), %rdx	#, tmp189
	leaq	stars(%rip), %rax	#, tmp190
	movss	(%rdx,%rax), %xmm1	# stars.vz, _25
# main.cpp:53:             stars.vz[i] += dz * stars.mass[j] * G_plus_dt * inv_d2;
	movq	32(%rsp), %rax	# j, tmp192
	addq	$288, %rax	#, tmp191
	leaq	0(,%rax,4), %rdx	#, tmp193
	leaq	stars(%rip), %rax	#, tmp194
	movss	(%rdx,%rax), %xmm0	# stars.mass, _26
# main.cpp:53:             stars.vz[i] += dz * stars.mass[j] * G_plus_dt * inv_d2;
	movaps	%xmm0, %xmm2	# _26, _26
	mulss	12(%rsp), %xmm2	# dz, _26
# main.cpp:53:             stars.vz[i] += dz * stars.mass[j] * G_plus_dt * inv_d2;
	movss	.LC3(%rip), %xmm0	#, tmp195
	mulss	%xmm2, %xmm0	# _27, _28
# main.cpp:53:             stars.vz[i] += dz * stars.mass[j] * G_plus_dt * inv_d2;
	mulss	20(%rsp), %xmm0	# inv_d2, _29
# main.cpp:53:             stars.vz[i] += dz * stars.mass[j] * G_plus_dt * inv_d2;
	addss	%xmm1, %xmm0	# _25, _30
	movq	24(%rsp), %rax	# i, tmp197
	addq	$240, %rax	#, tmp196
	leaq	0(,%rax,4), %rdx	#, tmp198
	leaq	stars(%rip), %rax	#, tmp199
	movss	%xmm0, (%rdx,%rax)	# _30, stars.vz
# main.cpp:44:         for (size_t j = 0; j < N; j++) {
	addq	$1, 32(%rsp)	#, j
	jmp	.L15	#
.L14:
# main.cpp:42:     for (size_t i = 0; i < N; i++) {
	addq	$1, 24(%rsp)	#, i
	jmp	.L16	#
.L13:
# main.cpp:57:     for (size_t i = 0; i < N; i++) {
	movq	$0, 40(%rsp)	#, i
.L18:
# main.cpp:57:     for (size_t i = 0; i < N; i++) {
	cmpq	$47, 40(%rsp)	#, i
	setbe	%al	#, _31
	testb	%al, %al	# _31
	je	.L19	#,
# main.cpp:58:         stars.px[i] += stars.vx[i] * dt;
	movq	40(%rsp), %rax	# i, tmp200
	leaq	0(,%rax,4), %rdx	#, tmp201
	leaq	stars(%rip), %rax	#, tmp202
	movss	(%rdx,%rax), %xmm1	# stars.px, _32
# main.cpp:58:         stars.px[i] += stars.vx[i] * dt;
	movq	40(%rsp), %rax	# i, tmp204
	addq	$144, %rax	#, tmp203
	leaq	0(,%rax,4), %rdx	#, tmp205
	leaq	stars(%rip), %rax	#, tmp206
	movss	(%rdx,%rax), %xmm2	# stars.vx, _33
# main.cpp:58:         stars.px[i] += stars.vx[i] * dt;
	movss	.LC4(%rip), %xmm0	#, tmp207
	mulss	%xmm2, %xmm0	# _33, _34
# main.cpp:58:         stars.px[i] += stars.vx[i] * dt;
	addss	%xmm1, %xmm0	# _32, _35
	movq	40(%rsp), %rax	# i, tmp208
	leaq	0(,%rax,4), %rdx	#, tmp209
	leaq	stars(%rip), %rax	#, tmp210
	movss	%xmm0, (%rdx,%rax)	# _35, stars.px
# main.cpp:59:         stars.py[i] += stars.vy[i] * dt;
	movq	40(%rsp), %rax	# i, tmp212
	addq	$48, %rax	#, tmp211
	leaq	0(,%rax,4), %rdx	#, tmp213
	leaq	stars(%rip), %rax	#, tmp214
	movss	(%rdx,%rax), %xmm1	# stars.py, _36
# main.cpp:59:         stars.py[i] += stars.vy[i] * dt;
	movq	40(%rsp), %rax	# i, tmp216
	addq	$192, %rax	#, tmp215
	leaq	0(,%rax,4), %rdx	#, tmp217
	leaq	stars(%rip), %rax	#, tmp218
	movss	(%rdx,%rax), %xmm2	# stars.vy, _37
# main.cpp:59:         stars.py[i] += stars.vy[i] * dt;
	movss	.LC4(%rip), %xmm0	#, tmp219
	mulss	%xmm2, %xmm0	# _37, _38
# main.cpp:59:         stars.py[i] += stars.vy[i] * dt;
	addss	%xmm1, %xmm0	# _36, _39
	movq	40(%rsp), %rax	# i, tmp221
	addq	$48, %rax	#, tmp220
	leaq	0(,%rax,4), %rdx	#, tmp222
	leaq	stars(%rip), %rax	#, tmp223
	movss	%xmm0, (%rdx,%rax)	# _39, stars.py
# main.cpp:60:         stars.pz[i] += stars.vz[i] * dt;
	movq	40(%rsp), %rax	# i, tmp225
	addq	$96, %rax	#, tmp224
	leaq	0(,%rax,4), %rdx	#, tmp226
	leaq	stars(%rip), %rax	#, tmp227
	movss	(%rdx,%rax), %xmm1	# stars.pz, _40
# main.cpp:60:         stars.pz[i] += stars.vz[i] * dt;
	movq	40(%rsp), %rax	# i, tmp229
	addq	$240, %rax	#, tmp228
	leaq	0(,%rax,4), %rdx	#, tmp230
	leaq	stars(%rip), %rax	#, tmp231
	movss	(%rdx,%rax), %xmm2	# stars.vz, _41
# main.cpp:60:         stars.pz[i] += stars.vz[i] * dt;
	movss	.LC4(%rip), %xmm0	#, tmp232
	mulss	%xmm2, %xmm0	# _41, _42
# main.cpp:60:         stars.pz[i] += stars.vz[i] * dt;
	addss	%xmm1, %xmm0	# _40, _43
	movq	40(%rsp), %rax	# i, tmp234
	addq	$96, %rax	#, tmp233
	leaq	0(,%rax,4), %rdx	#, tmp235
	leaq	stars(%rip), %rax	#, tmp236
	movss	%xmm0, (%rdx,%rax)	# _43, stars.pz
# main.cpp:57:     for (size_t i = 0; i < N; i++) {
	addq	$1, 40(%rsp)	#, i
	jmp	.L18	#
.L19:
# main.cpp:62: }
	nop	
	addq	$56, %rsp	#,
	.cfi_def_cfa_offset 8
	ret	
	.cfi_endproc
.LFE1371:
	.size	_Z4stepv, .-_Z4stepv
	.globl	_Z4calcv
	.type	_Z4calcv, @function
_Z4calcv:
.LFB1372:
	.cfi_startproc
	endbr64	
	subq	$56, %rsp	#,
	.cfi_def_cfa_offset 64
# main.cpp:65:     float energy = 0;
	pxor	%xmm0, %xmm0	# tmp115
	movss	%xmm0, (%rsp)	# tmp115, energy
# main.cpp:66:     float half_G = G * 0.5f;
	movss	.LC6(%rip), %xmm0	#, tmp116
	movss	%xmm0, 4(%rsp)	# tmp116, half_G
# main.cpp:67:     for (size_t i = 0; i < N; i++) {
	movq	$0, 32(%rsp)	#, i
.L24:
# main.cpp:67:     for (size_t i = 0; i < N; i++) {
	cmpq	$47, 32(%rsp)	#, i
	ja	.L21	#,
# main.cpp:68:         float v2 = stars.vx[i] * stars.vx[i] + stars.vy[i] * stars.vy[i] + stars.vz[i] * stars.vz[i];
	movq	32(%rsp), %rax	# i, tmp118
	addq	$144, %rax	#, tmp117
	leaq	0(,%rax,4), %rdx	#, tmp119
	leaq	stars(%rip), %rax	#, tmp120
	movss	(%rdx,%rax), %xmm1	# stars.vx, _1
# main.cpp:68:         float v2 = stars.vx[i] * stars.vx[i] + stars.vy[i] * stars.vy[i] + stars.vz[i] * stars.vz[i];
	movq	32(%rsp), %rax	# i, tmp122
	addq	$144, %rax	#, tmp121
	leaq	0(,%rax,4), %rdx	#, tmp123
	leaq	stars(%rip), %rax	#, tmp124
	movss	(%rdx,%rax), %xmm0	# stars.vx, _2
# main.cpp:68:         float v2 = stars.vx[i] * stars.vx[i] + stars.vy[i] * stars.vy[i] + stars.vz[i] * stars.vz[i];
	mulss	%xmm0, %xmm1	# _2, _3
# main.cpp:68:         float v2 = stars.vx[i] * stars.vx[i] + stars.vy[i] * stars.vy[i] + stars.vz[i] * stars.vz[i];
	movq	32(%rsp), %rax	# i, tmp126
	addq	$192, %rax	#, tmp125
	leaq	0(,%rax,4), %rdx	#, tmp127
	leaq	stars(%rip), %rax	#, tmp128
	movss	(%rdx,%rax), %xmm2	# stars.vy, _4
# main.cpp:68:         float v2 = stars.vx[i] * stars.vx[i] + stars.vy[i] * stars.vy[i] + stars.vz[i] * stars.vz[i];
	movq	32(%rsp), %rax	# i, tmp130
	addq	$192, %rax	#, tmp129
	leaq	0(,%rax,4), %rdx	#, tmp131
	leaq	stars(%rip), %rax	#, tmp132
	movss	(%rdx,%rax), %xmm0	# stars.vy, _5
# main.cpp:68:         float v2 = stars.vx[i] * stars.vx[i] + stars.vy[i] * stars.vy[i] + stars.vz[i] * stars.vz[i];
	mulss	%xmm2, %xmm0	# _4, _6
# main.cpp:68:         float v2 = stars.vx[i] * stars.vx[i] + stars.vy[i] * stars.vy[i] + stars.vz[i] * stars.vz[i];
	addss	%xmm0, %xmm1	# _6, _7
# main.cpp:68:         float v2 = stars.vx[i] * stars.vx[i] + stars.vy[i] * stars.vy[i] + stars.vz[i] * stars.vz[i];
	movq	32(%rsp), %rax	# i, tmp134
	addq	$240, %rax	#, tmp133
	leaq	0(,%rax,4), %rdx	#, tmp135
	leaq	stars(%rip), %rax	#, tmp136
	movss	(%rdx,%rax), %xmm2	# stars.vz, _8
# main.cpp:68:         float v2 = stars.vx[i] * stars.vx[i] + stars.vy[i] * stars.vy[i] + stars.vz[i] * stars.vz[i];
	movq	32(%rsp), %rax	# i, tmp138
	addq	$240, %rax	#, tmp137
	leaq	0(,%rax,4), %rdx	#, tmp139
	leaq	stars(%rip), %rax	#, tmp140
	movss	(%rdx,%rax), %xmm0	# stars.vz, _9
# main.cpp:68:         float v2 = stars.vx[i] * stars.vx[i] + stars.vy[i] * stars.vy[i] + stars.vz[i] * stars.vz[i];
	mulss	%xmm2, %xmm0	# _8, _10
# main.cpp:68:         float v2 = stars.vx[i] * stars.vx[i] + stars.vy[i] * stars.vy[i] + stars.vz[i] * stars.vz[i];
	addss	%xmm1, %xmm0	# _7, tmp141
	movss	%xmm0, 8(%rsp)	# tmp141, v2
# main.cpp:69:         energy += stars.mass[i] * v2 * 0.5f;
	movq	32(%rsp), %rax	# i, tmp143
	addq	$288, %rax	#, tmp142
	leaq	0(,%rax,4), %rdx	#, tmp144
	leaq	stars(%rip), %rax	#, tmp145
	movss	(%rdx,%rax), %xmm0	# stars.mass, _11
# main.cpp:69:         energy += stars.mass[i] * v2 * 0.5f;
	movaps	%xmm0, %xmm1	# _11, _11
	mulss	8(%rsp), %xmm1	# v2, _11
# main.cpp:69:         energy += stars.mass[i] * v2 * 0.5f;
	movss	.LC7(%rip), %xmm0	#, tmp146
	mulss	%xmm1, %xmm0	# _12, _13
# main.cpp:69:         energy += stars.mass[i] * v2 * 0.5f;
	movss	(%rsp), %xmm1	# energy, tmp148
	addss	%xmm1, %xmm0	# tmp148, tmp147
	movss	%xmm0, (%rsp)	# tmp147, energy
# main.cpp:71:         for (size_t j = 0; j < N; j++) {
	movq	$0, 40(%rsp)	#, j
.L23:
# main.cpp:71:         for (size_t j = 0; j < N; j++) {
	cmpq	$47, 40(%rsp)	#, j
	setbe	%al	#, _14
	testb	%al, %al	# _14
	je	.L22	#,
# main.cpp:72:             float dx = stars.px[j] - stars.px[i];
	movq	40(%rsp), %rax	# j, tmp149
	leaq	0(,%rax,4), %rdx	#, tmp150
	leaq	stars(%rip), %rax	#, tmp151
	movss	(%rdx,%rax), %xmm0	# stars.px, _15
# main.cpp:72:             float dx = stars.px[j] - stars.px[i];
	movq	32(%rsp), %rax	# i, tmp152
	leaq	0(,%rax,4), %rdx	#, tmp153
	leaq	stars(%rip), %rax	#, tmp154
	movss	(%rdx,%rax), %xmm1	# stars.px, _16
# main.cpp:72:             float dx = stars.px[j] - stars.px[i];
	subss	%xmm1, %xmm0	# _16, tmp155
	movss	%xmm0, 12(%rsp)	# tmp155, dx
# main.cpp:73:             float dy = stars.py[j] - stars.py[i];
	movq	40(%rsp), %rax	# j, tmp157
	addq	$48, %rax	#, tmp156
	leaq	0(,%rax,4), %rdx	#, tmp158
	leaq	stars(%rip), %rax	#, tmp159
	movss	(%rdx,%rax), %xmm0	# stars.py, _17
# main.cpp:73:             float dy = stars.py[j] - stars.py[i];
	movq	32(%rsp), %rax	# i, tmp161
	addq	$48, %rax	#, tmp160
	leaq	0(,%rax,4), %rdx	#, tmp162
	leaq	stars(%rip), %rax	#, tmp163
	movss	(%rdx,%rax), %xmm1	# stars.py, _18
# main.cpp:73:             float dy = stars.py[j] - stars.py[i];
	subss	%xmm1, %xmm0	# _18, tmp164
	movss	%xmm0, 16(%rsp)	# tmp164, dy
# main.cpp:74:             float dz = stars.pz[j] - stars.pz[i];
	movq	40(%rsp), %rax	# j, tmp166
	addq	$96, %rax	#, tmp165
	leaq	0(,%rax,4), %rdx	#, tmp167
	leaq	stars(%rip), %rax	#, tmp168
	movss	(%rdx,%rax), %xmm0	# stars.pz, _19
# main.cpp:74:             float dz = stars.pz[j] - stars.pz[i];
	movq	32(%rsp), %rax	# i, tmp170
	addq	$96, %rax	#, tmp169
	leaq	0(,%rax,4), %rdx	#, tmp171
	leaq	stars(%rip), %rax	#, tmp172
	movss	(%rdx,%rax), %xmm1	# stars.pz, _20
# main.cpp:74:             float dz = stars.pz[j] - stars.pz[i];
	subss	%xmm1, %xmm0	# _20, tmp173
	movss	%xmm0, 20(%rsp)	# tmp173, dz
# main.cpp:75:             float d2 = dx * dx + dy * dy + dz * dz + eps2;
	movss	12(%rsp), %xmm0	# dx, tmp174
	movaps	%xmm0, %xmm1	# tmp174, tmp174
	mulss	%xmm0, %xmm1	# tmp174, tmp174
# main.cpp:75:             float d2 = dx * dx + dy * dy + dz * dz + eps2;
	movss	16(%rsp), %xmm0	# dy, tmp175
	mulss	%xmm0, %xmm0	# tmp175, _22
# main.cpp:75:             float d2 = dx * dx + dy * dy + dz * dz + eps2;
	addss	%xmm0, %xmm1	# _22, _23
# main.cpp:75:             float d2 = dx * dx + dy * dy + dz * dz + eps2;
	movss	20(%rsp), %xmm0	# dz, tmp176
	mulss	%xmm0, %xmm0	# tmp176, _24
# main.cpp:75:             float d2 = dx * dx + dy * dy + dz * dz + eps2;
	addss	%xmm0, %xmm1	# _24, _25
# main.cpp:75:             float d2 = dx * dx + dy * dy + dz * dz + eps2;
	movss	.LC2(%rip), %xmm0	#, tmp178
	addss	%xmm1, %xmm0	# _25, tmp177
	movss	%xmm0, 24(%rsp)	# tmp177, d2
# main.cpp:76:             float inv_sqrt_d2 = 1.f / std::sqrt(d2);
	movl	24(%rsp), %eax	# d2, tmp179
	movd	%eax, %xmm0	# tmp179,
	call	_ZSt4sqrtf	#
# main.cpp:76:             float inv_sqrt_d2 = 1.f / std::sqrt(d2);
	movss	.LC1(%rip), %xmm1	#, tmp181
	divss	%xmm0, %xmm1	# _26, tmp181
	movaps	%xmm1, %xmm0	# tmp181, tmp180
	movss	%xmm0, 28(%rsp)	# tmp180, inv_sqrt_d2
# main.cpp:77:             energy -= stars.mass[j] * stars.mass[i] * half_G * inv_sqrt_d2;
	movq	40(%rsp), %rax	# j, tmp183
	addq	$288, %rax	#, tmp182
	leaq	0(,%rax,4), %rdx	#, tmp184
	leaq	stars(%rip), %rax	#, tmp185
	movss	(%rdx,%rax), %xmm1	# stars.mass, _27
# main.cpp:77:             energy -= stars.mass[j] * stars.mass[i] * half_G * inv_sqrt_d2;
	movq	32(%rsp), %rax	# i, tmp187
	addq	$288, %rax	#, tmp186
	leaq	0(,%rax,4), %rdx	#, tmp188
	leaq	stars(%rip), %rax	#, tmp189
	movss	(%rdx,%rax), %xmm0	# stars.mass, _28
# main.cpp:77:             energy -= stars.mass[j] * stars.mass[i] * half_G * inv_sqrt_d2;
	mulss	%xmm1, %xmm0	# _27, _29
# main.cpp:77:             energy -= stars.mass[j] * stars.mass[i] * half_G * inv_sqrt_d2;
	mulss	4(%rsp), %xmm0	# half_G, _30
# main.cpp:77:             energy -= stars.mass[j] * stars.mass[i] * half_G * inv_sqrt_d2;
	movaps	%xmm0, %xmm1	# _30, _30
	mulss	28(%rsp), %xmm1	# inv_sqrt_d2, _30
# main.cpp:77:             energy -= stars.mass[j] * stars.mass[i] * half_G * inv_sqrt_d2;
	movss	(%rsp), %xmm0	# energy, tmp191
	subss	%xmm1, %xmm0	# _31, tmp190
	movss	%xmm0, (%rsp)	# tmp190, energy
# main.cpp:71:         for (size_t j = 0; j < N; j++) {
	addq	$1, 40(%rsp)	#, j
	jmp	.L23	#
.L22:
# main.cpp:67:     for (size_t i = 0; i < N; i++) {
	addq	$1, 32(%rsp)	#, i
	jmp	.L24	#
.L21:
# main.cpp:80:     return energy;
	movss	(%rsp), %xmm0	# energy, _55
# main.cpp:81: }
	addq	$56, %rsp	#,
	.cfi_def_cfa_offset 8
	ret	
	.cfi_endproc
.LFE1372:
	.size	_Z4calcv, .-_Z4calcv
	.align 2
	.type	_ZZ4mainENKUlvE_clEv, @function
_ZZ4mainENKUlvE_clEv:
.LFB1375:
	.cfi_startproc
	subq	$40, %rsp	#,
	.cfi_def_cfa_offset 48
	movq	%rdi, 8(%rsp)	# __closure, __closure
# main.cpp:96:         for (int i = 0; i < 100000; i++)
	movl	$0, 28(%rsp)	#, i
.L28:
# main.cpp:96:         for (int i = 0; i < 100000; i++)
	cmpl	$99999, 28(%rsp)	#, i
	jg	.L29	#,
# main.cpp:97:             step();
	call	_Z4stepv	#
# main.cpp:96:         for (int i = 0; i < 100000; i++)
	addl	$1, 28(%rsp)	#, i
	jmp	.L28	#
.L29:
# main.cpp:98:     });
	nop	
	addq	$40, %rsp	#,
	.cfi_def_cfa_offset 8
	ret	
	.cfi_endproc
.LFE1375:
	.size	_ZZ4mainENKUlvE_clEv, .-_ZZ4mainENKUlvE_clEv
	.section	.rodata
.LC8:
	.string	"Initial energy: %f\n"
.LC9:
	.string	"Final energy: %f\n"
.LC10:
	.string	"Time elapsed: %ld ms\n"
	.text
	.globl	main
	.type	main, @function
main:
.LFB1374:
	.cfi_startproc
	endbr64	
	subq	$40, %rsp	#,
	.cfi_def_cfa_offset 48
# main.cpp:92: int main() {
	movq	%fs:40, %rax	# MEM[(<address-space-1> long unsigned int *)40B], tmp92
	movq	%rax, 24(%rsp)	# tmp92, D.28846
	xorl	%eax, %eax	# tmp92
# main.cpp:93:     init();
	call	_Z4initv	#
# main.cpp:94:     printf("Initial energy: %f\n", calc());
	call	_Z4calcv	#
# main.cpp:94:     printf("Initial energy: %f\n", calc());
	cvtss2sd	%xmm0, %xmm0	# _1, _2
	leaq	.LC8(%rip), %rdi	#,
	movl	$1, %eax	#,
	call	printf@PLT	#
# main.cpp:95:     auto dt = benchmark([&] {
	leaq	15(%rsp), %rax	#, tmp89
	movq	%rax, %rdi	# tmp89,
	call	_Z9benchmarkIZ4mainEUlvE_ElRKT_	#
	movq	%rax, 16(%rsp)	# _10, dt
# main.cpp:99:     printf("Final energy: %f\n", calc());
	call	_Z4calcv	#
# main.cpp:99:     printf("Final energy: %f\n", calc());
	cvtss2sd	%xmm0, %xmm0	# _3, _4
	leaq	.LC9(%rip), %rdi	#,
	movl	$1, %eax	#,
	call	printf@PLT	#
# main.cpp:100:     printf("Time elapsed: %ld ms\n", dt);
	movq	16(%rsp), %rax	# dt, tmp90
	movq	%rax, %rsi	# tmp90,
	leaq	.LC10(%rip), %rdi	#,
	movl	$0, %eax	#,
	call	printf@PLT	#
# main.cpp:101:     return 0;
	movl	$0, %eax	#, _16
# main.cpp:102: }
	movq	24(%rsp), %rdx	# D.28846, tmp93
	xorq	%fs:40, %rdx	# MEM[(<address-space-1> long unsigned int *)40B], tmp93
	je	.L32	#,
	call	__stack_chk_fail@PLT	#
.L32:
	addq	$40, %rsp	#,
	.cfi_def_cfa_offset 8
	ret	
	.cfi_endproc
.LFE1374:
	.size	main, .-main
	.section	.text._ZNSt6chrono13duration_castINS_8durationIlSt5ratioILl1ELl1000EEEElS2_ILl1ELl1000000000EEEENSt9enable_ifIXsrNS_13__is_durationIT_EE5valueES8_E4typeERKNS1_IT0_T1_EE,"axG",@progbits,_ZNSt6chrono13duration_castINS_8durationIlSt5ratioILl1ELl1000EEEElS2_ILl1ELl1000000000EEEENSt9enable_ifIXsrNS_13__is_durationIT_EE5valueES8_E4typeERKNS1_IT0_T1_EE,comdat
	.weak	_ZNSt6chrono13duration_castINS_8durationIlSt5ratioILl1ELl1000EEEElS2_ILl1ELl1000000000EEEENSt9enable_ifIXsrNS_13__is_durationIT_EE5valueES8_E4typeERKNS1_IT0_T1_EE
	.type	_ZNSt6chrono13duration_castINS_8durationIlSt5ratioILl1ELl1000EEEElS2_ILl1ELl1000000000EEEENSt9enable_ifIXsrNS_13__is_durationIT_EE5valueES8_E4typeERKNS1_IT0_T1_EE, @function
_ZNSt6chrono13duration_castINS_8durationIlSt5ratioILl1ELl1000EEEElS2_ILl1ELl1000000000EEEENSt9enable_ifIXsrNS_13__is_durationIT_EE5valueES8_E4typeERKNS1_IT0_T1_EE:
.LFB1407:
	.cfi_startproc
	endbr64	
	subq	$24, %rsp	#,
	.cfi_def_cfa_offset 32
	movq	%rdi, 8(%rsp)	# __d, __d
# /usr/include/c++/9/chrono:200: 	return __dc::__cast(__d);
	movq	8(%rsp), %rax	# __d, tmp84
	movq	%rax, %rdi	# tmp84,
	call	_ZNSt6chrono20__duration_cast_implINS_8durationIlSt5ratioILl1ELl1000EEEES2_ILl1ELl1000000EElLb1ELb0EE6__castIlS2_ILl1ELl1000000000EEEES4_RKNS1_IT_T0_EE	#
# /usr/include/c++/9/chrono:201:       }
	addq	$24, %rsp	#,
	.cfi_def_cfa_offset 8
	ret	
	.cfi_endproc
.LFE1407:
	.size	_ZNSt6chrono13duration_castINS_8durationIlSt5ratioILl1ELl1000EEEElS2_ILl1ELl1000000000EEEENSt9enable_ifIXsrNS_13__is_durationIT_EE5valueES8_E4typeERKNS1_IT0_T1_EE, .-_ZNSt6chrono13duration_castINS_8durationIlSt5ratioILl1ELl1000EEEElS2_ILl1ELl1000000000EEEENSt9enable_ifIXsrNS_13__is_durationIT_EE5valueES8_E4typeERKNS1_IT0_T1_EE
	.section	.text._ZNSt6chronomiINS_3_V212steady_clockENS_8durationIlSt5ratioILl1ELl1000000000EEEES6_EENSt11common_typeIJT0_T1_EE4typeERKNS_10time_pointIT_S8_EERKNSC_ISD_S9_EE,"axG",@progbits,_ZNSt6chronomiINS_3_V212steady_clockENS_8durationIlSt5ratioILl1ELl1000000000EEEES6_EENSt11common_typeIJT0_T1_EE4typeERKNS_10time_pointIT_S8_EERKNSC_ISD_S9_EE,comdat
	.weak	_ZNSt6chronomiINS_3_V212steady_clockENS_8durationIlSt5ratioILl1ELl1000000000EEEES6_EENSt11common_typeIJT0_T1_EE4typeERKNS_10time_pointIT_S8_EERKNSC_ISD_S9_EE
	.type	_ZNSt6chronomiINS_3_V212steady_clockENS_8durationIlSt5ratioILl1ELl1000000000EEEES6_EENSt11common_typeIJT0_T1_EE4typeERKNS_10time_pointIT_S8_EERKNSC_ISD_S9_EE, @function
_ZNSt6chronomiINS_3_V212steady_clockENS_8durationIlSt5ratioILl1ELl1000000000EEEES6_EENSt11common_typeIJT0_T1_EE4typeERKNS_10time_pointIT_S8_EERKNSC_ISD_S9_EE:
.LFB1408:
	.cfi_startproc
	endbr64	
	subq	$56, %rsp	#,
	.cfi_def_cfa_offset 64
	movq	%rdi, 8(%rsp)	# __lhs, __lhs
	movq	%rsi, (%rsp)	# __rhs, __rhs
# /usr/include/c++/9/chrono:762:       operator-(const time_point<_Clock, _Dur1>& __lhs,
	movq	%fs:40, %rax	# MEM[(<address-space-1> long unsigned int *)40B], tmp93
	movq	%rax, 40(%rsp)	# tmp93, D.28851
	xorl	%eax, %eax	# tmp93
# /usr/include/c++/9/chrono:764:       { return __lhs.time_since_epoch() - __rhs.time_since_epoch(); }
	movq	(%rsp), %rax	# __rhs, tmp84
	movq	%rax, %rdi	# tmp84,
	call	_ZNKSt6chrono10time_pointINS_3_V212steady_clockENS_8durationIlSt5ratioILl1ELl1000000000EEEEE16time_since_epochEv	#
	movq	%rax, 32(%rsp)	# tmp86, D.28631
	movq	8(%rsp), %rax	# __lhs, tmp87
	movq	%rax, %rdi	# tmp87,
	call	_ZNKSt6chrono10time_pointINS_3_V212steady_clockENS_8durationIlSt5ratioILl1ELl1000000000EEEEE16time_since_epochEv	#
	movq	%rax, 24(%rsp)	# tmp89, D.28630
	leaq	32(%rsp), %rdx	#, tmp90
	leaq	24(%rsp), %rax	#, tmp91
	movq	%rdx, %rsi	# tmp90,
	movq	%rax, %rdi	# tmp91,
	call	_ZNSt6chronomiIlSt5ratioILl1ELl1000000000EElS2_EENSt11common_typeIJNS_8durationIT_T0_EENS4_IT1_T2_EEEE4typeERKS7_RKSA_	#
# /usr/include/c++/9/chrono:764:       { return __lhs.time_since_epoch() - __rhs.time_since_epoch(); }
	movq	40(%rsp), %rcx	# D.28851, tmp94
	xorq	%fs:40, %rcx	# MEM[(<address-space-1> long unsigned int *)40B], tmp94
	je	.L37	#,
	call	__stack_chk_fail@PLT	#
.L37:
	addq	$56, %rsp	#,
	.cfi_def_cfa_offset 8
	ret	
	.cfi_endproc
.LFE1408:
	.size	_ZNSt6chronomiINS_3_V212steady_clockENS_8durationIlSt5ratioILl1ELl1000000000EEEES6_EENSt11common_typeIJT0_T1_EE4typeERKNS_10time_pointIT_S8_EERKNSC_ISD_S9_EE, .-_ZNSt6chronomiINS_3_V212steady_clockENS_8durationIlSt5ratioILl1ELl1000000000EEEES6_EENSt11common_typeIJT0_T1_EE4typeERKNS_10time_pointIT_S8_EERKNSC_ISD_S9_EE
	.section	.text._ZNSt6chronomiIlSt5ratioILl1ELl1000000000EElS2_EENSt11common_typeIJNS_8durationIT_T0_EENS4_IT1_T2_EEEE4typeERKS7_RKSA_,"axG",@progbits,_ZNSt6chronomiIlSt5ratioILl1ELl1000000000EElS2_EENSt11common_typeIJNS_8durationIT_T0_EENS4_IT1_T2_EEEE4typeERKS7_RKSA_,comdat
	.weak	_ZNSt6chronomiIlSt5ratioILl1ELl1000000000EElS2_EENSt11common_typeIJNS_8durationIT_T0_EENS4_IT1_T2_EEEE4typeERKS7_RKSA_
	.type	_ZNSt6chronomiIlSt5ratioILl1ELl1000000000EElS2_EENSt11common_typeIJNS_8durationIT_T0_EENS4_IT1_T2_EEEE4typeERKS7_RKSA_, @function
_ZNSt6chronomiIlSt5ratioILl1ELl1000000000EElS2_EENSt11common_typeIJNS_8durationIT_T0_EENS4_IT1_T2_EEEE4typeERKS7_RKSA_:
.LFB1409:
	.cfi_startproc
	endbr64	
	pushq	%rbx	#
	.cfi_def_cfa_offset 16
	.cfi_offset 3, -16
	subq	$64, %rsp	#,
	.cfi_def_cfa_offset 80
	movq	%rdi, 8(%rsp)	# __lhs, __lhs
	movq	%rsi, (%rsp)	# __rhs, __rhs
# /usr/include/c++/9/chrono:463:       operator-(const duration<_Rep1, _Period1>& __lhs,
	movq	%fs:40, %rax	# MEM[(<address-space-1> long unsigned int *)40B], tmp96
	movq	%rax, 56(%rsp)	# tmp96, D.28852
	xorl	%eax, %eax	# tmp96
# /usr/include/c++/9/chrono:469: 	return __cd(__cd(__lhs).count() - __cd(__rhs).count());
	movq	8(%rsp), %rax	# __lhs, tmp87
	movq	(%rax), %rax	# *__lhs_5(D), tmp88
	movq	%rax, 24(%rsp)	# tmp88, D.28647
	leaq	24(%rsp), %rax	#, tmp89
	movq	%rax, %rdi	# tmp89,
	call	_ZNKSt6chrono8durationIlSt5ratioILl1ELl1000000000EEE5countEv	#
	movq	%rax, %rbx	#, _1
	movq	(%rsp), %rax	# __rhs, tmp90
	movq	(%rax), %rax	# *__rhs_8(D), tmp91
	movq	%rax, 32(%rsp)	# tmp91, D.28649
	leaq	32(%rsp), %rax	#, tmp92
	movq	%rax, %rdi	# tmp92,
	call	_ZNKSt6chrono8durationIlSt5ratioILl1ELl1000000000EEE5countEv	#
	subq	%rax, %rbx	# _2, _1
	movq	%rbx, %rax	# _1, _3
	movq	%rax, 40(%rsp)	# _3, D.28650
# /usr/include/c++/9/chrono:469: 	return __cd(__cd(__lhs).count() - __cd(__rhs).count());
	leaq	40(%rsp), %rdx	#, tmp93
	leaq	48(%rsp), %rax	#, tmp94
	movq	%rdx, %rsi	# tmp93,
	movq	%rax, %rdi	# tmp94,
	call	_ZNSt6chrono8durationIlSt5ratioILl1ELl1000000000EEEC1IlvEERKT_	#
	movq	48(%rsp), %rax	# D.28651, D.28813
# /usr/include/c++/9/chrono:470:       }
	movq	56(%rsp), %rcx	# D.28852, tmp97
	xorq	%fs:40, %rcx	# MEM[(<address-space-1> long unsigned int *)40B], tmp97
	je	.L40	#,
	call	__stack_chk_fail@PLT	#
.L40:
	addq	$64, %rsp	#,
	.cfi_def_cfa_offset 16
	popq	%rbx	#
	.cfi_def_cfa_offset 8
	ret	
	.cfi_endproc
.LFE1409:
	.size	_ZNSt6chronomiIlSt5ratioILl1ELl1000000000EElS2_EENSt11common_typeIJNS_8durationIT_T0_EENS4_IT1_T2_EEEE4typeERKS7_RKSA_, .-_ZNSt6chronomiIlSt5ratioILl1ELl1000000000EElS2_EENSt11common_typeIJNS_8durationIT_T0_EENS4_IT1_T2_EEEE4typeERKS7_RKSA_
	.section	.text._ZNKSt6chrono10time_pointINS_3_V212steady_clockENS_8durationIlSt5ratioILl1ELl1000000000EEEEE16time_since_epochEv,"axG",@progbits,_ZNKSt6chrono10time_pointINS_3_V212steady_clockENS_8durationIlSt5ratioILl1ELl1000000000EEEEE16time_since_epochEv,comdat
	.align 2
	.weak	_ZNKSt6chrono10time_pointINS_3_V212steady_clockENS_8durationIlSt5ratioILl1ELl1000000000EEEEE16time_since_epochEv
	.type	_ZNKSt6chrono10time_pointINS_3_V212steady_clockENS_8durationIlSt5ratioILl1ELl1000000000EEEEE16time_since_epochEv, @function
_ZNKSt6chrono10time_pointINS_3_V212steady_clockENS_8durationIlSt5ratioILl1ELl1000000000EEEEE16time_since_epochEv:
.LFB1410:
	.cfi_startproc
	endbr64	
	movq	%rdi, -8(%rsp)	# this, this
# /usr/include/c++/9/chrono:650: 	{ return __d; }
	movq	-8(%rsp), %rax	# this, tmp84
	movq	(%rax), %rax	# this_2(D)->__d, D.28811
# /usr/include/c++/9/chrono:650: 	{ return __d; }
	ret	
	.cfi_endproc
.LFE1410:
	.size	_ZNKSt6chrono10time_pointINS_3_V212steady_clockENS_8durationIlSt5ratioILl1ELl1000000000EEEEE16time_since_epochEv, .-_ZNKSt6chrono10time_pointINS_3_V212steady_clockENS_8durationIlSt5ratioILl1ELl1000000000EEEEE16time_since_epochEv
	.text
	.type	_Z9benchmarkIZ4mainEUlvE_ElRKT_, @function
_Z9benchmarkIZ4mainEUlvE_ElRKT_:
.LFB1406:
	.cfi_startproc
	endbr64	
	subq	$72, %rsp	#,
	.cfi_def_cfa_offset 80
	movq	%rdi, 8(%rsp)	# func, func
# main.cpp:84: long benchmark(Func const &func) {
	movq	%fs:40, %rax	# MEM[(<address-space-1> long unsigned int *)40B], tmp98
	movq	%rax, 56(%rsp)	# tmp98, D.28853
	xorl	%eax, %eax	# tmp98
# main.cpp:85:     auto t0 = std::chrono::steady_clock::now();
	call	_ZNSt6chrono3_V212steady_clock3nowEv@PLT	#
	movq	%rax, 24(%rsp)	# tmp85, t0
# main.cpp:86:     func();
	movq	8(%rsp), %rax	# func, tmp86
	movq	%rax, %rdi	# tmp86,
	call	_ZZ4mainENKUlvE_clEv	#
# main.cpp:87:     auto t1 = std::chrono::steady_clock::now();
	call	_ZNSt6chrono3_V212steady_clock3nowEv@PLT	#
	movq	%rax, 32(%rsp)	# tmp88, t1
# main.cpp:88:     auto dt = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0);
	leaq	24(%rsp), %rdx	#, tmp89
	leaq	32(%rsp), %rax	#, tmp90
	movq	%rdx, %rsi	# tmp89,
	movq	%rax, %rdi	# tmp90,
	call	_ZNSt6chronomiINS_3_V212steady_clockENS_8durationIlSt5ratioILl1ELl1000000000EEEES6_EENSt11common_typeIJT0_T1_EE4typeERKNS_10time_pointIT_S8_EERKNSC_ISD_S9_EE	#
	movq	%rax, 48(%rsp)	# tmp92, D.28567
# main.cpp:88:     auto dt = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0);
	leaq	48(%rsp), %rax	#, tmp93
	movq	%rax, %rdi	# tmp93,
	call	_ZNSt6chrono13duration_castINS_8durationIlSt5ratioILl1ELl1000EEEElS2_ILl1ELl1000000000EEEENSt9enable_ifIXsrNS_13__is_durationIT_EE5valueES8_E4typeERKNS1_IT0_T1_EE	#
	movq	%rax, 40(%rsp)	# tmp95, dt
# main.cpp:89:     return dt.count();
	leaq	40(%rsp), %rax	#, tmp96
	movq	%rax, %rdi	# tmp96,
	call	_ZNKSt6chrono8durationIlSt5ratioILl1ELl1000EEE5countEv	#
# main.cpp:90: }
	movq	56(%rsp), %rcx	# D.28853, tmp99
	xorq	%fs:40, %rcx	# MEM[(<address-space-1> long unsigned int *)40B], tmp99
	je	.L45	#,
	call	__stack_chk_fail@PLT	#
.L45:
	addq	$72, %rsp	#,
	.cfi_def_cfa_offset 8
	ret	
	.cfi_endproc
.LFE1406:
	.size	_Z9benchmarkIZ4mainEUlvE_ElRKT_, .-_Z9benchmarkIZ4mainEUlvE_ElRKT_
	.section	.text._ZNSt6chrono20__duration_cast_implINS_8durationIlSt5ratioILl1ELl1000EEEES2_ILl1ELl1000000EElLb1ELb0EE6__castIlS2_ILl1ELl1000000000EEEES4_RKNS1_IT_T0_EE,"axG",@progbits,_ZNSt6chrono20__duration_cast_implINS_8durationIlSt5ratioILl1ELl1000EEEES2_ILl1ELl1000000EElLb1ELb0EE6__castIlS2_ILl1ELl1000000000EEEES4_RKNS1_IT_T0_EE,comdat
	.weak	_ZNSt6chrono20__duration_cast_implINS_8durationIlSt5ratioILl1ELl1000EEEES2_ILl1ELl1000000EElLb1ELb0EE6__castIlS2_ILl1ELl1000000000EEEES4_RKNS1_IT_T0_EE
	.type	_ZNSt6chrono20__duration_cast_implINS_8durationIlSt5ratioILl1ELl1000EEEES2_ILl1ELl1000000EElLb1ELb0EE6__castIlS2_ILl1ELl1000000000EEEES4_RKNS1_IT_T0_EE, @function
_ZNSt6chrono20__duration_cast_implINS_8durationIlSt5ratioILl1ELl1000EEEES2_ILl1ELl1000000EElLb1ELb0EE6__castIlS2_ILl1ELl1000000000EEEES4_RKNS1_IT_T0_EE:
.LFB1417:
	.cfi_startproc
	endbr64	
	subq	$56, %rsp	#,
	.cfi_def_cfa_offset 64
	movq	%rdi, 8(%rsp)	# __d, __d
# /usr/include/c++/9/chrono:149: 	  __cast(const duration<_Rep, _Period>& __d)
	movq	%fs:40, %rax	# MEM[(<address-space-1> long unsigned int *)40B], tmp94
	movq	%rax, 40(%rsp)	# tmp94, D.28854
	xorl	%eax, %eax	# tmp94
# /usr/include/c++/9/chrono:153: 	      static_cast<_CR>(__d.count()) / static_cast<_CR>(_CF::den)));
	movq	8(%rsp), %rax	# __d, tmp86
	movq	%rax, %rdi	# tmp86,
	call	_ZNKSt6chrono8durationIlSt5ratioILl1ELl1000000000EEE5countEv	#
	movq	%rax, %rcx	#, _1
# /usr/include/c++/9/chrono:153: 	      static_cast<_CR>(__d.count()) / static_cast<_CR>(_CF::den)));
	movabsq	$4835703278458516699, %rdx	#, tmp88
	movq	%rcx, %rax	# _1, tmp95
	imulq	%rdx	# tmp88
	sarq	$18, %rdx	#, tmp89
	movq	%rcx, %rax	# _1, _1
	sarq	$63, %rax	#, _1
	subq	%rax, %rdx	# tmp90, tmp89
	movq	%rdx, %rax	# tmp89, _2
	movq	%rax, 24(%rsp)	# _2, D.28748
# /usr/include/c++/9/chrono:153: 	      static_cast<_CR>(__d.count()) / static_cast<_CR>(_CF::den)));
	leaq	24(%rsp), %rdx	#, tmp91
	leaq	32(%rsp), %rax	#, tmp92
	movq	%rdx, %rsi	# tmp91,
	movq	%rax, %rdi	# tmp92,
	call	_ZNSt6chrono8durationIlSt5ratioILl1ELl1000EEEC1IlvEERKT_	#
	movq	32(%rsp), %rax	# D.28749, D.28823
# /usr/include/c++/9/chrono:154: 	  }
	movq	40(%rsp), %rsi	# D.28854, tmp96
	xorq	%fs:40, %rsi	# MEM[(<address-space-1> long unsigned int *)40B], tmp96
	je	.L48	#,
	call	__stack_chk_fail@PLT	#
.L48:
	addq	$56, %rsp	#,
	.cfi_def_cfa_offset 8
	ret	
	.cfi_endproc
.LFE1417:
	.size	_ZNSt6chrono20__duration_cast_implINS_8durationIlSt5ratioILl1ELl1000EEEES2_ILl1ELl1000000EElLb1ELb0EE6__castIlS2_ILl1ELl1000000000EEEES4_RKNS1_IT_T0_EE, .-_ZNSt6chrono20__duration_cast_implINS_8durationIlSt5ratioILl1ELl1000EEEES2_ILl1ELl1000000EElLb1ELb0EE6__castIlS2_ILl1ELl1000000000EEEES4_RKNS1_IT_T0_EE
	.section	.text._ZNKSt6chrono8durationIlSt5ratioILl1ELl1000EEE5countEv,"axG",@progbits,_ZNKSt6chrono8durationIlSt5ratioILl1ELl1000EEE5countEv,comdat
	.align 2
	.weak	_ZNKSt6chrono8durationIlSt5ratioILl1ELl1000EEE5countEv
	.type	_ZNKSt6chrono8durationIlSt5ratioILl1ELl1000EEE5countEv, @function
_ZNKSt6chrono8durationIlSt5ratioILl1ELl1000EEE5countEv:
.LFB1418:
	.cfi_startproc
	endbr64	
	movq	%rdi, -8(%rsp)	# this, this
# /usr/include/c++/9/chrono:347: 	{ return __r; }
	movq	-8(%rsp), %rax	# this, tmp84
	movq	(%rax), %rax	# this_2(D)->__r, _3
# /usr/include/c++/9/chrono:347: 	{ return __r; }
	ret	
	.cfi_endproc
.LFE1418:
	.size	_ZNKSt6chrono8durationIlSt5ratioILl1ELl1000EEE5countEv, .-_ZNKSt6chrono8durationIlSt5ratioILl1ELl1000EEE5countEv
	.section	.text._ZNSt6chrono8durationIlSt5ratioILl1ELl1000EEEC2IlvEERKT_,"axG",@progbits,_ZNSt6chrono8durationIlSt5ratioILl1ELl1000EEEC5IlvEERKT_,comdat
	.align 2
	.weak	_ZNSt6chrono8durationIlSt5ratioILl1ELl1000EEEC2IlvEERKT_
	.type	_ZNSt6chrono8durationIlSt5ratioILl1ELl1000EEEC2IlvEERKT_, @function
_ZNSt6chrono8durationIlSt5ratioILl1ELl1000EEEC2IlvEERKT_:
.LFB1421:
	.cfi_startproc
	endbr64	
	movq	%rdi, -8(%rsp)	# this, this
	movq	%rsi, -16(%rsp)	# __rep, __rep
# /usr/include/c++/9/chrono:332: 	  : __r(static_cast<rep>(__rep)) { }
	movq	-16(%rsp), %rax	# __rep, tmp83
	movq	(%rax), %rdx	# *__rep_5(D), _1
	movq	-8(%rsp), %rax	# this, tmp84
	movq	%rdx, (%rax)	# _1, this_3(D)->__r
# /usr/include/c++/9/chrono:332: 	  : __r(static_cast<rep>(__rep)) { }
	nop	
	ret	
	.cfi_endproc
.LFE1421:
	.size	_ZNSt6chrono8durationIlSt5ratioILl1ELl1000EEEC2IlvEERKT_, .-_ZNSt6chrono8durationIlSt5ratioILl1ELl1000EEEC2IlvEERKT_
	.weak	_ZNSt6chrono8durationIlSt5ratioILl1ELl1000EEEC1IlvEERKT_
	.set	_ZNSt6chrono8durationIlSt5ratioILl1ELl1000EEEC1IlvEERKT_,_ZNSt6chrono8durationIlSt5ratioILl1ELl1000EEEC2IlvEERKT_
	.section	.rodata
	.align 4
.LC0:
	.long	1325400064
	.align 4
.LC1:
	.long	1065353216
	.align 4
.LC2:
	.long	897988542
	.align 4
.LC3:
	.long	925353389
	.align 4
.LC4:
	.long	1008981770
	.align 4
.LC6:
	.long	973279855
	.align 4
.LC7:
	.long	1056964608
	.ident	"GCC: (Ubuntu 9.4.0-1ubuntu1~20.04.1) 9.4.0"
	.section	.note.GNU-stack,"",@progbits
	.section	.note.gnu.property,"a"
	.align 8
	.long	 1f - 0f
	.long	 4f - 1f
	.long	 5
0:
	.string	 "GNU"
1:
	.align 8
	.long	 0xc0000002
	.long	 3f - 2f
2:
	.long	 0x3
3:
	.align 8
4:
