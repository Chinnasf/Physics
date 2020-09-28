! DATE: MAR. 2019
! AUTHORS: KARINA CHINAS F. & DANIEL HERNANDEZ M.
! CODE SIMULATES A QUANTUM TELEPORTATION OF 4 QUBITS
! THE SYSTEM IS BASED ON A CARBON CRYSTAL STRUCTURE WITH C13 PARTICLES (AS QUBITS) 


program Four_Qubit
    implicit none
    !VARIABLE DECLARATION
    integer, parameter :: out_unit=20
    real, parameter:: PI=3.1415922653589793238
    character(len=20):: file_name

    integer:: n,i
    real :: ti, tf, t, tf1, tf2, tf3, tf4, tf5, hpp, pp
    real :: h, U, theta, w, phi, w0, w1, w2, w3, w4, J1 , J2, J3
    real :: Energy, E01, E02, E03, E04, E05, E06, E07, E08, E09, E10, E11, E12, E13, E14, E15, E16

    real :: B01_02, B01_03, B01_05, B01_09, B02_04, B02_06, B02_10, B03_04
    real :: B03_07, B03_11, B04_08, B04_12, B05_06, B05_07, B05_13, B06_08
    real :: B06_14, B07_08, B07_15, B08_16, B09_10, B09_11, B09_13, B10_12
    real :: B10_14, B11_12, B11_15, B12_16, B13_14, B13_15, B14_16, B15_16


    real :: DX01, DX02, DX03, DX04, DX05, DX06, DX07, DX08, DX09, DX10, DX11, DX12, DX13, DX14, DX15, DX16
    real :: DY01, DY02, DY03, DY04, DY05, DY06, DY07, DY08, DY09, DY10, DY11, DY12, DY13, DY14, DY15, DY16

    real :: X01, X02, X03, X04, X05, X06, X07, X08, X09, X10, X11, X12, X13, X14, X15, X16
    real :: Y01, Y02, Y03, Y04, Y05, Y06, Y07, Y08, Y09, Y10, Y11, Y12, Y13, Y14, Y15, Y16

    real :: CX01, CX02, CX03, CX04, CX05, CX06, CX07, CX08, CX09, CX10, CX11, CX12, CX13, CX14, CX15, CX16
    real :: CY01, CY02, CY03, CY04, CY05, CY06, CY07, CY08, CY09, CY10, CY11, CY12, CY13, CY14, CY15, CY16

    real :: D01, D02, D03, D04, D05, D06, D07, D08, D09, D10, D11, D12, D13, D14, D15, D16, TP

    real :: K01a, K01b, K01c, K01d
    real :: L01a, L01b, L01c, L01d

    real :: K02a, K02b, K02c, K02d
    real :: L02a, L02b, L02c, L02d

    real :: K03a, K03b, K03c, K03d
    real :: L03a, L03b, L03c, L03d

    real :: K04a, K04b, K04c, K04d
    real :: L04a, L04b, L04c, L04d

    real :: K05a, K05b, K05c, K05d
    real :: L05a, L05b, L05c, L05d

    real :: K06a, K06b, K06c, K06d
    real :: L06a, L06b, L06c, L06d

    real :: K07a, K07b, K07c, K07d
    real :: L07a, L07b, L07c, L07d

    real :: K08a, K08b, K08c, K08d
    real :: L08a, L08b, L08c, L08d

    real :: K09a, K09b, K09c, K09d
    real :: L09a, L09b, L09c, L09d

    real :: K10a, K10b, K10c, K10d
    real :: L10a, L10b, L10c, L10d

    real :: K11a, K11b, K11c, K11d
    real :: L11a, L11b, L11c, L11d

    real :: K12a, K12b, K12c, K12d
    real :: L12a, L12b, L12c, L12d

    real :: K13a, K13b, K13c, K13d
    real :: L13a, L13b, L13c, L13d

    real :: K14a, K14b, K14c, K14d
    real :: L14a, L14b, L14c, L14d

    real :: K15a, K15b, K15c, K15d
    real :: L15a, L15b, L15c, L15d

    real :: K16a, K16b, K16c, K16d
    real :: L16a, L16b, L16c, L16d

!Global CONSTANTS FOR THE WHOLE SCRIPT

    common U, w, w0, phi, w1, w2, w3, w4, J1, J2, J3, E01, E02, E03, E04, &
    E05, E06, E07, E08, E09, E10, E11, E12, E13, E14, E15, E16

    w0  = 1.0/PI  !SCALED FREQUENCY FOR THE TIME TAU
    U   = 1.0/w0  !RABI FRECUENCY
    w1  = 100./w0 !FIRST-QUBIT LARMOR FREQUENCY
    w2  = 150./w0 !SECOND-QUBIT LARMOR FREQUENCY
    w3  = 250./w0 !THIRD-QUBIT LARMOR FREQUENCY
    w4  = 300./w0 !FOURTH-QUBIT LARMOR FREQUENCY
    J1  = 5.0/w0  !INTERACTION FREQUENCY OF FIRST NEIGHBORS
    J2  = J1/10.0 !INTERACTION FREQUENCY OF SECOND NEIGHBORS
    J3  = J2/10.0 !INTERACTION FREQUENCY OF THIRD NEIGHBORS

    phi = 0.0   ! MAGNETIC FIELD PHASE
	ti = 0.0	! INITIAL TIME - SECONDS: always zero
    
    pp = PI/U
    hpp = PI/(2*U)

    tf1 = hpp
    tf2 = tf1 + pp
    tf3 = tf2 + pp
    tf4 = tf3 + hpp
    tf5 = tf4 + hpp

    tf  = tf5 + hpp	! FINAL TIME

    h  = 0.001! STEP SIZE
    n  = int((tf-ti)/h) ! TOTAL NUMBER OF ITERATIONS
    i  = 0

    ! ENERGIES OF THE SYSTEM FOR EACH STATE

    E01 = Energy(0.0, 0.0, 0.0, 0.0)
    E02 = Energy(0.0, 0.0, 0.0, 1.0)
    E03 = Energy(0.0, 0.0, 1.0, 0.0)
    E04 = Energy(0.0, 0.0, 1.0, 1.0)
    E05 = Energy(0.0, 1.0, 0.0, 0.0)
    E06 = Energy(0.0, 1.0, 0.0, 1.0)
    E07 = Energy(0.0, 1.0, 1.0, 0.0)
    E08 = Energy(0.0, 1.0, 1.0, 1.0)
    E09 = Energy(1.0, 0.0, 0.0, 0.0)
    E10 = Energy(1.0, 0.0, 0.0, 1.0)
    E11 = Energy(1.0, 0.0, 1.0, 0.0)
    E12 = Energy(1.0, 0.0, 1.0, 1.0)
    E13 = Energy(1.0, 1.0, 0.0, 0.0)
    E14 = Energy(1.0, 1.0, 0.0, 1.0)
    E15 = Energy(1.0, 1.0, 1.0, 0.0)
    E16 = Energy(1.0, 1.0, 1.0, 1.0)

    w = E07-E05 !MAGNETIC FIELD FRECUENCY

    !INITIAL CONDITIONS OF THE QUANTUM SYSYEM

    t = ti

    !INITIAL CONDICIONS OF THE COEFICIENTS FOR THE INTERACTION SCHEME

    X01 = sqrt(1.0/2.0)
    Y01 = 0.0

    X02 = 0.0
    Y02 = 0.0
    X03 = 0.0
    Y03 = 0.0
    X04 = 0.0
    Y04 = 0.0
    X05 = 0.0
    Y05 = 0.0
    X06 = 0.0
    Y06 = 0.0
    X07 = 0.0
    Y07 = 0.0
    X08 = 0.0
    Y08 = 0.0
    
    X09 = 0.0
    Y09 = sqrt(1.0/2.0)
    
    X10 = 0.0
    Y10 = 0.0
    X11 = 0.0
    Y11 = 0.0
    X12 = 0.0
    Y12 = 0.0
    X13 = 0.0
    Y13 = 0.0
    X14 = 0.0
    Y14 = 0.0
    X15 = 0.0
    Y15 = 0.0
    X16 = 0.0
    Y16 = 0.0

    ! INITIAL CONDITIONS OF THE COEFICIENTS FOR THE WAVE FUNCTION
    CX01 = X01
    CY01 = Y01
    CX02 = X02
    CY02 = Y02
    CX03 = X03
    CY03 = Y03
    CX04 = X04
    CY04 = Y04
    CX05 = X05
    CY05 = Y05
    CX06 = X06
    CY06 = Y06
    CX07 = X07
    CY07 = Y07
    CX08 = X08
    CY08 = Y08
    CX09 = X09
    CY09 = Y09
    CX10 = X10
    CY10 = Y10
    CX11 = X11
    CY11 = Y11
    CX12 = X12
    CY12 = Y12
    CX13 = X13
    CY13 = Y13
    CX14 = X14
    CY14 = Y14
    CX15 = X15
    CY15 = Y15
    CX16 = X16
    CY16 = Y16

    ! PROBABILITY FOR EACH STATE AND TOTAL PROBABILITY

    D01 = (CX01**2) + (CY01**2)
    D02 = (CX02**2) + (CY02**2)
    D03 = (CX03**2) + (CY03**2)
    D04 = (CX04**2) + (CY04**2)
    D05 = (CX05**2) + (CY05**2)
    D06 = (CX06**2) + (CY06**2)
    D07 = (CX07**2) + (CY07**2)
    D08 = (CX08**2) + (CY08**2)
    D09 = (CX09**2) + (CY09**2)
    D10 = (CX10**2) + (CY10**2)
    D11 = (CX11**2) + (CY11**2)
    D12 = (CX12**2) + (CY12**2)
    D13 = (CX13**2) + (CY13**2)
    D14 = (CX14**2) + (CY14**2)
    D15 = (CX15**2) + (CY15**2)
    D16 = (CX16**2) + (CY16**2)

    TP = D01 + D02 + D03 + D04 + D05 + D06 + D07 + D08 + D09 + D10 + D11 + D12 + D13 + D14 + D15 + D16

    ! GENERATE A .TXT FILE TO RECORD ALL VALUES

    file_name = "ReQT4Q_R1.txt"		!MEANING OF THE NAME: "RESULTS 4-qubit Quantum Teleportation RUN number 1"
    open (unit=out_unit, file = file_name , action="write", status="replace")

    write(out_unit,*) t, X01, X02,   X03,  X04,  X05,  X06,  X07,  X08,  X09,  X10,  X11,  X12,  X13,  X14,  X15,  X16, &
                         Y01, Y02,   Y03,  Y04,  Y05,  Y06,  Y07,  Y08,  Y09,  Y10,  Y11,  Y12,  Y13,  Y14,  Y15,  Y16, &
                         D01 , D02 , D03 , D04 , D05 , D06 , D07 , D08 , D09 , D10 , D11 , D12 , D13 , D14 , D15 , D16 , TP
    do i = 1,n

        t = t + h

        if (t.le.tf1) then
            w = E02 - E01

        end if

        if (t.ge.tf1.and.t.le.tf2) then
            w = E06 - E02

        end if

        if (t.ge.tf2.and.t.le.tf3) then 
            w = E13 - E09

        end if    

        if (t.ge.tf3.and.t.le.tf4) then
            w = E14 - E06

        end if

        if (t.ge.tf4.and.t.le.tf5) then
            w = E10 - E02

        end if

        if (t.gt.tf5.and.t.le.tf) then
            w = 0

        end if 

        if (t.gt.tf) goto 500

        K01a = h*dX01(t, X02, Y02, X03, Y03, X05, Y05, X09, Y09)
        K02a = h*dX02(t, X01, Y01, X04, Y04, X06, Y06, X10, Y10)
        K03a = h*dX03(t, X01, Y01, X04, Y04, X07, Y07, X11, Y11)
        K04a = h*dX04(t, X02, Y02, X03, Y03, X08, Y08, X12, Y12)
        K05a = h*dX05(t, X01, Y01, X06, Y06, X07, Y07, X13, Y13)
        K06a = h*dX06(t, X02, Y02, X05, Y05, X08, Y08, X14, Y14)
        K07a = h*dX07(t, X03, Y03, X05, Y05, X08, Y08, X15, Y15)
        K08a = h*dX08(t, X04, Y04, X06, Y06, X07, Y07, X16, Y16)
        K09a = h*dX09(t, X01, Y01, X10, Y10, X11, Y11, X13, Y13)
        K10a = h*dX10(t, X02, Y02, X09, Y09, X12, Y12, X14, Y14)
        K11a = h*dX11(t, X03, Y03, X09, Y09, X12, Y12, X15, Y15)
        K12a = h*dX12(t, X04, Y04, X10, Y10, X11, Y11, X16, Y16)
        K13a = h*dX13(t, X05, Y05, X09, Y09, X14, Y14, X15, Y15)
        K14a = h*dX14(t, X06, Y06, X10, Y10, X13, Y13, X16, Y16)
        K15a = h*dX15(t, X07, Y07, X11, Y11, X13, Y13, X16, Y16)
        K16a = h*dX16(t, X08, Y08, X12, Y12, X14, Y14, X15, Y15)

        L01a = h*dY01(t, X02, Y02, X03, Y03, X05, Y05, X09, Y09)
        L02a = h*dY02(t, X01, Y01, X04, Y04, X06, Y06, X10, Y10)
        L03a = h*dY03(t, X01, Y01, X04, Y04, X07, Y07, X11, Y11)
        L04a = h*dY04(t, X02, Y02, X03, Y03, X08, Y08, X12, Y12)
        L05a = h*dY05(t, X01, Y01, X06, Y06, X07, Y07, X13, Y13)
        L06a = h*dY06(t, X02, Y02, X05, Y05, X08, Y08, X14, Y14)
        L07a = h*dY07(t, X03, Y03, X05, Y05, X08, Y08, X15, Y15)
        L08a = h*dY08(t, X04, Y04, X06, Y06, X07, Y07, X16, Y16)
        L09a = h*dY09(t, X01, Y01, X10, Y10, X11, Y11, X13, Y13)
        L10a = h*dY10(t, X02, Y02, X09, Y09, X12, Y12, X14, Y14)
        L11a = h*dY11(t, X03, Y03, X09, Y09, X12, Y12, X15, Y15)
        L12a = h*dY12(t, X04, Y04, X10, Y10, X11, Y11, X16, Y16)
        L13a = h*dY13(t, X05, Y05, X09, Y09, X14, Y14, X15, Y15)
        L14a = h*dY14(t, X06, Y06, X10, Y10, X13, Y13, X16, Y16)
        L15a = h*dY15(t, X07, Y07, X11, Y11, X13, Y13, X16, Y16)
        L16a = h*dY16(t, X08, Y08, X12, Y12, X14, Y14, X15, Y15)


        K01b = h*dX01(t + h*0.5, X02 + K02a*0.5, Y02 + L02a*0.5, X03 + K03a*0.5, &
               Y03 + L03a*0.5, X05 + K05a*0.5, Y05 + L05a*0.5, X09 + K09a*0.5, Y09 + L09a*0.5)

        K02b = h*dX02(t + h*0.5, X01 + K01a*0.5, Y01 + L01a*0.5, X04 + K04a*0.5, &
               Y04 + L04a*0.5, X06 + K06a*0.5, Y06 + L06a*0.5, X10 + K10a*0.5, Y10 + L10a*0.5)

        K03b = h*dX03(t + h*0.5, X01 + K01a*0.5, Y01 + L01a*0.5, X04 + K04a*0.5, &
               Y04 + L04a*0.5, X07 + K07a*0.5, Y07 + L07a*0.5, X11 + K11a*0.5, Y11 + L11a*0.5)

        K04b = h*dX04(t + h*0.5, X02 + K02a*0.5, Y02 + L02a*0.5, X03 + K03a*0.5, &
               Y03 + L03a*0.5, X08 + K08a*0.5, Y08 + L08a*0.5, X12 + K12a*0.5, Y12 + L12a*0.5)

        K05b = h*dX05(t + h*0.5, X01 + K01a*0.5, Y01 + L01a*0.5, X06 + K06a*0.5, &
               Y06 + L06a*0.5, X07 + K07a*0.5, Y07 + L07a*0.5, X13 + K13a*0.5, Y13 + L13a*0.5)

        K06b = h*dX06(t + h*0.5, X02 + K02a*0.5, Y02 + L02a*0.5, X05 + K05a*0.5, &
               Y05 + L05a*0.5, X08 + K08a*0.5, Y08 + L08a*0.5, X14 + K14a*0.5, Y14 + L14a*0.5)

        K07b = h*dX07(t + h*0.5, X03 + K03a*0.5, Y03 + L03a*0.5, X05 + K05a*0.5, &
               Y05 + L05a*0.5, X08 + K08a*0.5, Y08 + L08a*0.5, X15 + K15a*0.5, Y15 + L15a*0.5)

        K08b = h*dX08(t + h*0.5, X04 + K04a*0.5, Y04 + L04a*0.5, X06 + K06a*0.5, &
               Y06 + L06a*0.5, X07 + K07a*0.5, Y07 + L07a*0.5, X16 + K16a*0.5, Y16 + L16a*0.5)

        K09b = h*dX09(t + h*0.5, X01 + K01a*0.5, Y01 + L01a*0.5, X10 + K10a*0.5, &
               Y10 + L10a*0.5, X11 + K11a*0.5, Y11 + L11a*0.5, X13 + K13a*0.5, Y13 + L13a*0.5)

        K10b = h*dX10(t + h*0.5, X02 + K02a*0.5, Y02 + L02a*0.5, X09 + K09a*0.5, &
               Y09 + L09a*0.5, X12 + K12a*0.5, Y12 + L12a*0.5, X14 + K14a*0.5, Y14 + L14a*0.5)

        K11b = h*dX11(t + h*0.5, X03 + K03a*0.5, Y03 + L03a*0.5, X09 + K09a*0.5, &
               Y09 + L09a*0.5, X12 + K12a*0.5, Y12 + L12a*0.5, X15 + K15a*0.5, Y15 + L15a*0.5)

        K12b = h*dX12(t + h*0.5, X04 + K04a*0.5, Y04 + L04a*0.5, X10 + K10a*0.5, &
               Y10 + L10a*0.5, X11 + K11a*0.5, Y11 + L11a*0.5, X16 + K16a*0.5, Y16 + L16a*0.5)

        K13b = h*dX13(t + h*0.5, X05 + K05a*0.5, Y05 + L05a*0.5, X09 + K09a*0.5, &
               Y09 + L09a*0.5, X14 + K14a*0.5, Y14 + L14a*0.5, X15 + K15a*0.5, Y15 + L15a*0.5)

        K14b = h*dX14(t + h*0.5, X06 + K06a*0.5, Y06 + L06a*0.5, X10 + K10a*0.5, &
               Y10 + L10a*0.5, X13 + K13a*0.5, Y13 + L13a*0.5, X16 + K16a*0.5, Y16 + L16a*0.5)

        K15b = h*dX15(t + h*0.5, X07 + K07a*0.5, Y07 + L07a*0.5, X11 + K11a*0.5, &
               Y11 + L11a*0.5, X13 + K13a*0.5, Y13 + L13a*0.5, X16 + K16a*0.5, Y16 + L16a*0.5)

        K16b = h*dX16(t + h*0.5, X08 + K08a*0.5, Y08 + L08a*0.5, X12 + K12a*0.5, &
               Y12 + L12a*0.5, X14 + K14a*0.5, Y14 + L14a*0.5, X15 + K15a*0.5, Y15 + L15a*0.5)


        L01b = h*dY01(t + h*0.5, X02 + K02a*0.5, Y02 + L02a*0.5, X03 + K03a*0.5, & 
               Y03 + L03a*0.5, X05 + K05a*0.5, Y05 + L05a*0.5, X09 + K09a*0.5, Y09 + L09a*0.5)

        L02b = h*dY02(t + h*0.5, X01 + K01a*0.5, Y01 + L01a*0.5, X04 + K04a*0.5, & 
               Y04 + L04a*0.5, X06 + K06a*0.5, Y06 + L06a*0.5, X10 + K10a*0.5, Y10 + L10a*0.5)

        L03b = h*dY03(t + h*0.5, X01 + K01a*0.5, Y01 + L01a*0.5, X04 + K04a*0.5, & 
               Y04 + L04a*0.5, X07 + K07a*0.5, Y07 + L07a*0.5, X11 + K11a*0.5, Y11 + L11a*0.5)

        L04b = h*dY04(t + h*0.5, X02 + K02a*0.5, Y02 + L02a*0.5, X03 + K03a*0.5, & 
               Y03 + L03a*0.5, X08 + K08a*0.5, Y08 + L08a*0.5, X12 + K12a*0.5, Y12 + L12a*0.5)

        L05b = h*dY05(t + h*0.5, X01 + K01a*0.5, Y01 + L01a*0.5, X06 + K06a*0.5, & 
               Y06 + L06a*0.5, X07 + K07a*0.5, Y07 + L07a*0.5, X13 + K13a*0.5, Y13 + L13a*0.5)

        L06b = h*dY06(t + h*0.5, X02 + K02a*0.5, Y02 + L02a*0.5, X05 + K05a*0.5, & 
               Y05 + L05a*0.5, X08 + K08a*0.5, Y08 + L08a*0.5, X14 + K14a*0.5, Y14 + L14a*0.5)

        L07b = h*dY07(t + h*0.5, X03 + K03a*0.5, Y03 + L03a*0.5, X05 + K05a*0.5, & 
               Y05 + L05a*0.5, X08 + K08a*0.5, Y08 + L08a*0.5, X15 + K15a*0.5, Y15 + L15a*0.5)

        L08b = h*dY08(t + h*0.5, X04 + K04a*0.5, Y04 + L04a*0.5, X06 + K06a*0.5, & 
               Y06 + L06a*0.5, X07 + K07a*0.5, Y07 + L07a*0.5, X16 + K16a*0.5, Y16 + L16a*0.5)

        L09b = h*dY09(t + h*0.5, X01 + K01a*0.5, Y01 + L01a*0.5, X10 + K10a*0.5, & 
               Y10 + L10a*0.5, X11 + K11a*0.5, Y11 + L11a*0.5, X13 + K13a*0.5, Y13 + L13a*0.5)

        L10b = h*dY10(t + h*0.5, X02 + K02a*0.5, Y02 + L02a*0.5, X09 + K09a*0.5, & 
               Y09 + L09a*0.5, X12 + K12a*0.5, Y12 + L12a*0.5, X14 + K14a*0.5, Y14 + L14a*0.5)

        L11b = h*dY11(t + h*0.5, X03 + K03a*0.5, Y03 + L03a*0.5, X09 + K09a*0.5, & 
               Y09 + L09a*0.5, X12 + K12a*0.5, Y12 + L12a*0.5, X15 + K15a*0.5, Y15 + L15a*0.5)

        L12b = h*dY12(t + h*0.5, X04 + K04a*0.5, Y04 + L04a*0.5, X10 + K10a*0.5, & 
               Y10 + L10a*0.5, X11 + K11a*0.5, Y11 + L11a*0.5, X16 + K16a*0.5, Y16 + L16a*0.5)

        L13b = h*dY13(t + h*0.5, X05 + K05a*0.5, Y05 + L05a*0.5, X09 + K09a*0.5, & 
               Y09 + L09a*0.5, X14 + K14a*0.5, Y14 + L14a*0.5, X15 + K15a*0.5, Y15 + L15a*0.5)

        L14b = h*dY14(t + h*0.5, X06 + K06a*0.5, Y06 + L06a*0.5, X10 + K10a*0.5, & 
               Y10 + L10a*0.5, X13 + K13a*0.5, Y13 + L13a*0.5, X16 + K16a*0.5, Y16 + L16a*0.5)

        L15b = h*dY15(t + h*0.5, X07 + K07a*0.5, Y07 + L07a*0.5, X11 + K11a*0.5, & 
               Y11 + L11a*0.5, X13 + K13a*0.5, Y13 + L13a*0.5, X16 + K16a*0.5, Y16 + L16a*0.5)

        L16b = h*dY16(t + h*0.5, X08 + K08a*0.5, Y08 + L08a*0.5, X12 + K12a*0.5, & 
               Y12 + L12a*0.5, X14 + K14a*0.5, Y14 + L14a*0.5, X15 + K15a*0.5, Y15 + L15a*0.5)



        K01c = h*dX01(t + h*0.5, X02 + K02b*0.5, Y02 + L02b*0.5, X03 + & 
               K03b*0.5, Y03 + L03b*0.5, X05 + K05b*0.5, Y05 + L05b*0.5, X09 + K09b*0.5, Y09 + L09b*0.5)

        K02c = h*dX02(t + h*0.5, X01 + K01b*0.5, Y01 + L01b*0.5, X04 + & 
               K04b*0.5, Y04 + L04b*0.5, X06 + K06b*0.5, Y06 + L06b*0.5, X10 + K10b*0.5, Y10 + L10b*0.5)

        K03c = h*dX03(t + h*0.5, X01 + K01b*0.5, Y01 + L01b*0.5, X04 + & 
               K04b*0.5, Y04 + L04b*0.5, X07 + K07b*0.5, Y07 + L07b*0.5, X11 + K11b*0.5, Y11 + L11b*0.5)

        K04c = h*dX04(t + h*0.5, X02 + K02b*0.5, Y02 + L02b*0.5, X03 + & 
               K03b*0.5, Y03 + L03b*0.5, X08 + K08b*0.5, Y08 + L08b*0.5, X12 + K12b*0.5, Y12 + L12b*0.5)

        K05c = h*dX05(t + h*0.5, X01 + K01b*0.5, Y01 + L01b*0.5, X06 + & 
               K06b*0.5, Y06 + L06b*0.5, X07 + K07b*0.5, Y07 + L07b*0.5, X13 + K13b*0.5, Y13 + L13b*0.5)

        K06c = h*dX06(t + h*0.5, X02 + K02b*0.5, Y02 + L02b*0.5, X05 + & 
               K05b*0.5, Y05 + L05b*0.5, X08 + K08b*0.5, Y08 + L08b*0.5, X14 + K14b*0.5, Y14 + L14b*0.5)

        K07c = h*dX07(t + h*0.5, X03 + K03b*0.5, Y03 + L03b*0.5, X05 + & 
               K05b*0.5, Y05 + L05b*0.5, X08 + K08b*0.5, Y08 + L08b*0.5, X15 + K15b*0.5, Y15 + L15b*0.5)

        K08c = h*dX08(t + h*0.5, X04 + K04b*0.5, Y04 + L04b*0.5, X06 + & 
               K06b*0.5, Y06 + L06b*0.5, X07 + K07b*0.5, Y07 + L07b*0.5, X16 + K16b*0.5, Y16 + L16b*0.5)

        K09c = h*dX09(t + h*0.5, X01 + K01b*0.5, Y01 + L01b*0.5, X10 + & 
               K10b*0.5, Y10 + L10b*0.5, X11 + K11b*0.5, Y11 + L11b*0.5, X13 + K13b*0.5, Y13 + L13b*0.5)

        K10c = h*dX10(t + h*0.5, X02 + K02b*0.5, Y02 + L02b*0.5, X09 + & 
               K09b*0.5, Y09 + L09b*0.5, X12 + K12b*0.5, Y12 + L12b*0.5, X14 + K14b*0.5, Y14 + L14b*0.5)

        K11c = h*dX11(t + h*0.5, X03 + K03b*0.5, Y03 + L03b*0.5, X09 + & 
               K09b*0.5, Y09 + L09b*0.5, X12 + K12b*0.5, Y12 + L12b*0.5, X15 + K15b*0.5, Y15 + L15b*0.5)

        K12c = h*dX12(t + h*0.5, X04 + K04b*0.5, Y04 + L04b*0.5, X10 + & 
               K10b*0.5, Y10 + L10b*0.5, X11 + K11b*0.5, Y11 + L11b*0.5, X16 + K16b*0.5, Y16 + L16b*0.5)

        K13c = h*dX13(t + h*0.5, X05 + K05b*0.5, Y05 + L05b*0.5, X09 + & 
               K09b*0.5, Y09 + L09b*0.5, X14 + K14b*0.5, Y14 + L14b*0.5, X15 + K15b*0.5, Y15 + L15b*0.5)

        K14c = h*dX14(t + h*0.5, X06 + K06b*0.5, Y06 + L06b*0.5, X10 + & 
               K10b*0.5, Y10 + L10b*0.5, X13 + K13b*0.5, Y13 + L13b*0.5, X16 + K16b*0.5, Y16 + L16b*0.5)

        K15c = h*dX15(t + h*0.5, X07 + K07b*0.5, Y07 + L07b*0.5, X11 + & 
               K11b*0.5, Y11 + L11b*0.5, X13 + K13b*0.5, Y13 + L13b*0.5, X16 + K16b*0.5, Y16 + L16b*0.5)

        K16c = h*dX16(t + h*0.5, X08 + K08b*0.5, Y08 + L08b*0.5, X12 + & 
               K12b*0.5, Y12 + L12b*0.5, X14 + K14b*0.5, Y14 + L14b*0.5, X15 + K15b*0.5, Y15 + L15b*0.5)


        L01c = h*dY01(t + h*0.5, X02 + K02b*0.5, Y02 + L02b*0.5, X03 + & 
               K03b*0.5, Y03 + L03b*0.5, X05 + K05b*0.5, Y05 + L05b*0.5, X09 + K09b*0.5, Y09 + L09b*0.5)

        L02c = h*dY02(t + h*0.5, X01 + K01b*0.5, Y01 + L01b*0.5, X04 + & 
               K04b*0.5, Y04 + L04b*0.5, X06 + K06b*0.5, Y06 + L06b*0.5, X10 + K10b*0.5, Y10 + L10b*0.5)

        L03c = h*dY03(t + h*0.5, X01 + K01b*0.5, Y01 + L01b*0.5, X04 + & 
               K04b*0.5, Y04 + L04b*0.5, X07 + K07b*0.5, Y07 + L07b*0.5, X11 + K11b*0.5, Y11 + L11b*0.5)

        L04c = h*dY04(t + h*0.5, X02 + K02b*0.5, Y02 + L02b*0.5, X03 + & 
               K03b*0.5, Y03 + L03b*0.5, X08 + K08b*0.5, Y08 + L08b*0.5, X12 + K12b*0.5, Y12 + L12b*0.5)

        L05c = h*dY05(t + h*0.5, X01 + K01b*0.5, Y01 + L01b*0.5, X06 + & 
               K06b*0.5, Y06 + L06b*0.5, X07 + K07b*0.5, Y07 + L07b*0.5, X13 + K13b*0.5, Y13 + L13b*0.5)

        L06c = h*dY06(t + h*0.5, X02 + K02b*0.5, Y02 + L02b*0.5, X05 + & 
               K05b*0.5, Y05 + L05b*0.5, X08 + K08b*0.5, Y08 + L08b*0.5, X14 + K14b*0.5, Y14 + L14b*0.5)

        L07c = h*dY07(t + h*0.5, X03 + K03b*0.5, Y03 + L03b*0.5, X05 + & 
               K05b*0.5, Y05 + L05b*0.5, X08 + K08b*0.5, Y08 + L08b*0.5, X15 + K15b*0.5, Y15 + L15b*0.5)

        L08c = h*dY08(t + h*0.5, X04 + K04b*0.5, Y04 + L04b*0.5, X06 + & 
               K06b*0.5, Y06 + L06b*0.5, X07 + K07b*0.5, Y07 + L07b*0.5, X16 + K16b*0.5, Y16 + L16b*0.5)

        L09c = h*dY09(t + h*0.5, X01 + K01b*0.5, Y01 + L01b*0.5, X10 + & 
               K10b*0.5, Y10 + L10b*0.5, X11 + K11b*0.5, Y11 + L11b*0.5, X13 + K13b*0.5, Y13 + L13b*0.5)

        L10c = h*dY10(t + h*0.5, X02 + K02b*0.5, Y02 + L02b*0.5, X09 + & 
               K09b*0.5, Y09 + L09b*0.5, X12 + K12b*0.5, Y12 + L12b*0.5, X14 + K14b*0.5, Y14 + L14b*0.5)

        L11c = h*dY11(t + h*0.5, X03 + K03b*0.5, Y03 + L03b*0.5, X09 + & 
               K09b*0.5, Y09 + L09b*0.5, X12 + K12b*0.5, Y12 + L12b*0.5, X15 + K15b*0.5, Y15 + L15b*0.5)

        L12c = h*dY12(t + h*0.5, X04 + K04b*0.5, Y04 + L04b*0.5, X10 + & 
               K10b*0.5, Y10 + L10b*0.5, X11 + K11b*0.5, Y11 + L11b*0.5, X16 + K16b*0.5, Y16 + L16b*0.5)

        L13c = h*dY13(t + h*0.5, X05 + K05b*0.5, Y05 + L05b*0.5, X09 + & 
               K09b*0.5, Y09 + L09b*0.5, X14 + K14b*0.5, Y14 + L14b*0.5, X15 + K15b*0.5, Y15 + L15b*0.5)

        L14c = h*dY14(t + h*0.5, X06 + K06b*0.5, Y06 + L06b*0.5, X10 + & 
               K10b*0.5, Y10 + L10b*0.5, X13 + K13b*0.5, Y13 + L13b*0.5, X16 + K16b*0.5, Y16 + L16b*0.5)

        L15c = h*dY15(t + h*0.5, X07 + K07b*0.5, Y07 + L07b*0.5, X11 + & 
               K11b*0.5, Y11 + L11b*0.5, X13 + K13b*0.5, Y13 + L13b*0.5, X16 + K16b*0.5, Y16 + L16b*0.5)

        L16c = h*dY16(t + h*0.5, X08 + K08b*0.5, Y08 + L08b*0.5, X12 + & 
               K12b*0.5, Y12 + L12b*0.5, X14 + K14b*0.5, Y14 + L14b*0.5, X15 + K15b*0.5, Y15 + L15b*0.5)



        K01d = h*dX01(t + h, X02 + K02c, Y02 + L02c, X03 + K03c, Y03 + L03c, X05 + K05c, Y05 + L05c, X09 + K09c, Y09 + L09c)
        K02d = h*dX02(t + h, X01 + K01c, Y01 + L01c, X04 + K04c, Y04 + L04c, X06 + K06c, Y06 + L06c, X10 + K10c, Y10 + L10c)
        K03d = h*dX03(t + h, X01 + K01c, Y01 + L01c, X04 + K04c, Y04 + L04c, X07 + K07c, Y07 + L07c, X11 + K11c, Y11 + L11c)
        K04d = h*dX04(t + h, X02 + K02c, Y02 + L02c, X03 + K03c, Y03 + L03c, X08 + K08c, Y08 + L08c, X12 + K12c, Y12 + L12c)
        K05d = h*dX05(t + h, X01 + K01c, Y01 + L01c, X06 + K06c, Y06 + L06c, X07 + K07c, Y07 + L07c, X13 + K13c, Y13 + L13c)
        K06d = h*dX06(t + h, X02 + K02c, Y02 + L02c, X05 + K05c, Y05 + L05c, X08 + K08c, Y08 + L08c, X14 + K14c, Y14 + L14c)
        K07d = h*dX07(t + h, X03 + K03c, Y03 + L03c, X05 + K05c, Y05 + L05c, X08 + K08c, Y08 + L08c, X15 + K15c, Y15 + L15c)
        K08d = h*dX08(t + h, X04 + K04c, Y04 + L04c, X06 + K06c, Y06 + L06c, X07 + K07c, Y07 + L07c, X16 + K16c, Y16 + L16c)
        K09d = h*dX09(t + h, X01 + K01c, Y01 + L01c, X10 + K10c, Y10 + L10c, X11 + K11c, Y11 + L11c, X13 + K13c, Y13 + L13c)
        K10d = h*dX10(t + h, X02 + K02c, Y02 + L02c, X09 + K09c, Y09 + L09c, X12 + K12c, Y12 + L12c, X14 + K14c, Y14 + L14c)
        K11d = h*dX11(t + h, X03 + K03c, Y03 + L03c, X09 + K09c, Y09 + L09c, X12 + K12c, Y12 + L12c, X15 + K15c, Y15 + L15c)
        K12d = h*dX12(t + h, X04 + K04c, Y04 + L04c, X10 + K10c, Y10 + L10c, X11 + K11c, Y11 + L11c, X16 + K16c, Y16 + L16c)
        K13d = h*dX13(t + h, X05 + K05c, Y05 + L05c, X09 + K09c, Y09 + L09c, X14 + K14c, Y14 + L14c, X15 + K15c, Y15 + L15c)
        K14d = h*dX14(t + h, X06 + K06c, Y06 + L06c, X10 + K10c, Y10 + L10c, X13 + K13c, Y13 + L13c, X16 + K16c, Y16 + L16c)
        K15d = h*dX15(t + h, X07 + K07c, Y07 + L07c, X11 + K11c, Y11 + L11c, X13 + K13c, Y13 + L13c, X16 + K16c, Y16 + L16c)
        K16d = h*dX16(t + h, X08 + K08c, Y08 + L08c, X12 + K12c, Y12 + L12c, X14 + K14c, Y14 + L14c, X15 + K15c, Y15 + L15c)

        L01d = h*dY01(t + h, X02 + K02c, Y02 + L02c, X03 + K03c, Y03 + L03c, X05 + K05c, Y05 + L05c, X09 + K09c, Y09 + L09c)
        L02d = h*dY02(t + h, X01 + K01c, Y01 + L01c, X04 + K04c, Y04 + L04c, X06 + K06c, Y06 + L06c, X10 + K10c, Y10 + L10c)
        L03d = h*dY03(t + h, X01 + K01c, Y01 + L01c, X04 + K04c, Y04 + L04c, X07 + K07c, Y07 + L07c, X11 + K11c, Y11 + L11c)
        L04d = h*dY04(t + h, X02 + K02c, Y02 + L02c, X03 + K03c, Y03 + L03c, X08 + K08c, Y08 + L08c, X12 + K12c, Y12 + L12c)
        L05d = h*dY05(t + h, X01 + K01c, Y01 + L01c, X06 + K06c, Y06 + L06c, X07 + K07c, Y07 + L07c, X13 + K13c, Y13 + L13c)
        L06d = h*dY06(t + h, X02 + K02c, Y02 + L02c, X05 + K05c, Y05 + L05c, X08 + K08c, Y08 + L08c, X14 + K14c, Y14 + L14c)
        L07d = h*dY07(t + h, X03 + K03c, Y03 + L03c, X05 + K05c, Y05 + L05c, X08 + K08c, Y08 + L08c, X15 + K15c, Y15 + L15c)
        L08d = h*dY08(t + h, X04 + K04c, Y04 + L04c, X06 + K06c, Y06 + L06c, X07 + K07c, Y07 + L07c, X16 + K16c, Y16 + L16c)
        L09d = h*dY09(t + h, X01 + K01c, Y01 + L01c, X10 + K10c, Y10 + L10c, X11 + K11c, Y11 + L11c, X13 + K13c, Y13 + L13c)
        L10d = h*dY10(t + h, X02 + K02c, Y02 + L02c, X09 + K09c, Y09 + L09c, X12 + K12c, Y12 + L12c, X14 + K14c, Y14 + L14c)
        L11d = h*dY11(t + h, X03 + K03c, Y03 + L03c, X09 + K09c, Y09 + L09c, X12 + K12c, Y12 + L12c, X15 + K15c, Y15 + L15c)
        L12d = h*dY12(t + h, X04 + K04c, Y04 + L04c, X10 + K10c, Y10 + L10c, X11 + K11c, Y11 + L11c, X16 + K16c, Y16 + L16c)
        L13d = h*dY13(t + h, X05 + K05c, Y05 + L05c, X09 + K09c, Y09 + L09c, X14 + K14c, Y14 + L14c, X15 + K15c, Y15 + L15c)
        L14d = h*dY14(t + h, X06 + K06c, Y06 + L06c, X10 + K10c, Y10 + L10c, X13 + K13c, Y13 + L13c, X16 + K16c, Y16 + L16c)
        L15d = h*dY15(t + h, X07 + K07c, Y07 + L07c, X11 + K11c, Y11 + L11c, X13 + K13c, Y13 + L13c, X16 + K16c, Y16 + L16c)
        L16d = h*dY16(t + h, X08 + K08c, Y08 + L08c, X12 + K12c, Y12 + L12c, X14 + K14c, Y14 + L14c, X15 + K15c, Y15 + L15c)

        ! COEFICIENTS OF THE INTERACTION SCHEME (WITH LOW FREQUENCIES)

        X01 = X01 + (1.0/6.0)*( K01a + 2.0*K01b + 2.0*K01c + K01d )
        Y01 = Y01 + (1.0/6.0)*( L01a + 2.0*L01b + 2.0*L01c + L01d )

        X02 = X02 + (1.0/6.0)*( K02a + 2.0*K02b + 2.0*K02c + K02d )
        Y02 = Y02 + (1.0/6.0)*( L02a + 2.0*L02b + 2.0*L02c + L02d )

        X03 = X03 + (1.0/6.0)*( K03a + 2.0*K03b + 2.0*K03c + K03d )
        Y03 = Y03 + (1.0/6.0)*( L03a + 2.0*L03b + 2.0*L03c + L03d )

        X04 = X04 + (1.0/6.0)*( K04a + 2.0*K04b + 2.0*K04c + K04d )
        Y04 = Y04 + (1.0/6.0)*( L04a + 2.0*L04b + 2.0*L04c + L04d )

        X05 = X05 + (1.0/6.0)*( K05a + 2.0*K05b + 2.0*K05c + K05d )
        Y05 = Y05 + (1.0/6.0)*( L05a + 2.0*L05b + 2.0*L05c + L05d )

        X06 = X06 + (1.0/6.0)*( K06a + 2.0*K06b + 2.0*K06c + K06d )
        Y06 = Y06 + (1.0/6.0)*( L06a + 2.0*L06b + 2.0*L06c + L06d )

        X07 = X07 + (1.0/6.0)*( K07a + 2.0*K07b + 2.0*K07c + K07d )
        Y07 = Y07 + (1.0/6.0)*( L07a + 2.0*L07b + 2.0*L07c + L07d )

        X08 = X08 + (1.0/6.0)*( K08a + 2.0*K08b + 2.0*K08c + K08d )
        Y08 = Y08 + (1.0/6.0)*( L08a + 2.0*L08b + 2.0*L08c + L08d )

        X09 = X09 + (1.0/6.0)*( K09a + 2.0*K09b + 2.0*K09c + K09d )
        Y09 = Y09 + (1.0/6.0)*( L09a + 2.0*L09b + 2.0*L09c + L09d )

        X10 = X10 + (1.0/6.0)*( K10a + 2.0*K10b + 2.0*K10c + K10d )
        Y10 = Y10 + (1.0/6.0)*( L10a + 2.0*L10b + 2.0*L10c + L10d )

        X11 = X11 + (1.0/6.0)*( K11a + 2.0*K11b + 2.0*K11c + K11d )
        Y11 = Y11 + (1.0/6.0)*( L11a + 2.0*L11b + 2.0*L11c + L11d )

        X12 = X12 + (1.0/6.0)*( K12a + 2.0*K12b + 2.0*K12c + K12d )
        Y12 = Y12 + (1.0/6.0)*( L12a + 2.0*L12b + 2.0*L12c + L12d )

        X13 = X13 + (1.0/6.0)*( K13a + 2.0*K13b + 2.0*K13c + K13d )
        Y13 = Y13 + (1.0/6.0)*( L13a + 2.0*L13b + 2.0*L13c + L13d )

        X14 = X14 + (1.0/6.0)*( K14a + 2.0*K14b + 2.0*K14c + K14d )
        Y14 = Y14 + (1.0/6.0)*( L14a + 2.0*L14b + 2.0*L14c + L14d )

        X15 = X15 + (1.0/6.0)*( K15a + 2.0*K15b + 2.0*K15c + K15d )
        Y15 = Y15 + (1.0/6.0)*( L15a + 2.0*L15b + 2.0*L15c + L15d )

        X16 = X16 + (1.0/6.0)*( K16a + 2.0*K16b + 2.0*K16c + K16d )
        Y16 = Y16 + (1.0/6.0)*( L16a + 2.0*L16b + 2.0*L16c + L16d )

        ! COEFICIENTS OF THE WAVE FUNCTION (WITH HIGH FREQUENCIES)

        CX01 = X01*cos(E01*t) + (-Y01*sin(E01*t))
        CY01 = X01*sin(E01*t) +   Y01*cos(E01*t)

        CX02 = X02*cos(E02*t) + (-Y02*sin(E02*t))
        CY02 = X02*sin(E02*t) +   Y02*cos(E02*t)

        CX03 = X03*cos(E03*t) + (-Y03*sin(E03*t))
        CY03 = X03*sin(E03*t) +   Y03*cos(E03*t)

        CX04 = X04*cos(E04*t) + (-Y04*sin(E04*t))
        CY04 = X04*sin(E04*t) +   Y04*cos(E04*t)

        CX05 = X05*cos(E05*t) + (-Y05*sin(E05*t))
        CY05 = X05*sin(E05*t) +   Y05*cos(E05*t)

        CX06 = X06*cos(E06*t) + (-Y06*sin(E06*t))
        CY06 = X06*sin(E06*t) +   Y06*cos(E06*t)

        CX07 = X07*cos(E07*t) + (-Y07*sin(E07*t))
        CY07 = X07*sin(E07*t) +   Y07*cos(E07*t)

        CX08 = X08*cos(E08*t) + (-Y08*sin(E08*t))
        CY08 = X08*sin(E08*t) +   Y08*cos(E08*t)

        CX09 = X09*cos(E09*t) + (-Y09*sin(E09*t))
        CY09 = X09*sin(E09*t) +   Y09*cos(E09*t)

        CX10 = X10*cos(E10*t) + (-Y10*sin(E10*t))
        CY10 = X10*sin(E10*t) +   Y10*cos(E10*t)

        CX11 = X11*cos(E11*t) + (-Y11*sin(E11*t))
        CY11 = X11*sin(E11*t) +   Y11*cos(E11*t)

        CX12 = X12*cos(E12*t) + (-Y12*sin(E12*t))
        CY12 = X12*sin(E12*t) +   Y12*cos(E12*t)

        CX13 = X13*cos(E13*t) + (-Y13*sin(E13*t))
        CY13 = X13*sin(E13*t) +   Y13*cos(E13*t)

        CX14 = X14*cos(E14*t) + (-Y14*sin(E14*t))
        CY14 = X14*sin(E14*t) +   Y14*cos(E14*t)

        CX15 = X15*cos(E15*t) + (-Y15*sin(E15*t))
        CY15 = X15*sin(E15*t) +   Y15*cos(E15*t)

        CX16 = X16*cos(E16*t) + (-Y16*sin(E16*t))
        CY16 = X16*sin(E16*t) +   Y16*cos(E16*t)

        ! PROBABILITY OF EACH STATE AND TOTAL PROBABILITY

        D01 = (CX01**2) + (CY01**2)
        D02 = (CX02**2) + (CY02**2)
        D03 = (CX03**2) + (CY03**2)
        D04 = (CX04**2) + (CY04**2)
        D05 = (CX05**2) + (CY05**2)
        D06 = (CX06**2) + (CY06**2)
        D07 = (CX07**2) + (CY07**2)
        D08 = (CX08**2) + (CY08**2)
        D09 = (CX09**2) + (CY09**2)
        D10 = (CX10**2) + (CY10**2)
        D11 = (CX11**2) + (CY11**2)
        D12 = (CX12**2) + (CY12**2)
        D13 = (CX13**2) + (CY13**2)
        D14 = (CX14**2) + (CY14**2)
        D15 = (CX15**2) + (CY15**2)
        D16 = (CX16**2) + (CY16**2)

        TP = D01 + D02 + D03 + D04 + D05 + D06 + D07 + D08 + D09 + D10 + D11 + D12 + D13 + D14 + D15 + D16

        write(out_unit,*) t, X01, X02,   X03,  X04,  X05,  X06,  X07,  X08,  X09,  X10,  X11,  X12,  X13,  X14,  X15,  X16, &
                             Y01, Y02,   Y03,  Y04,  Y05,  Y06,  Y07,  Y08,  Y09,  Y10,  Y11,  Y12,  Y13,  Y14,  Y15,  Y16, &
                             D01 , D02 , D03 , D04 , D05 , D06 , D07 , D08 , D09 , D10 , D11 , D12 , D13 , D14 , D15 , D16 , TP
    end do
    500 continue
    close (out_unit)
end program Four_Qubit

! SCALED ENERGY OVER w0 AND h-BAR

real function Energy(a4, a3, a2, a1)
    implicit none

    real:: U, w, w0, phi,  w1, w2, w3, w4, J1, J2, J3, a1, a2, a3, a4, E01, E02, E03, E04, & 
    E05, E06, E07, E08, E09, E10, E11, E12, E13, E14, E15, E16

    common U, w, w0, phi, w1, w2, w3, w4, J1, J2, J3, E01, E02, E03, E04, &
    E05, E06, E07, E08, E09, E10, E11, E12, E13, E14, E15, E16

    Energy  = (-0.5/w0)*( w1*((-1.0)**a1) + w2*((-1.0)**a2) + w3*((-1)**a3) + w4*((-1)**a4) + &
          0.5*J1*(((-1.0)**(a1+a2)) + ((-1.0)**(a2+a3)) + ((-1.0)**(a3+a4)) ) + &
          0.5*J2*(((-1.0)**(a1+a3)) + ((-1.0)**(a2+a4)) ) + &
          0.5*J3*(((-1.0)**(a1+a4))  ) )
end function

! Differential equations

real function dX01(t, X02, Y02, X03, Y03, X05, Y05, X09, Y09)


    real :: U, theta, w, w0, phi,  w1, w2, w3, w4, J1, J2, J3
    real :: E01, E02, E03, E04, E05, E06, E07, E08, E09, E10, E11, E12, E13, E14, E15, E16
    real :: t, X02, Y02, X03, Y03, X05, Y05, X09, Y09

    real :: B01_02, B01_03, B01_05, B01_09, B02_04, B02_06, B02_10, B03_04
    real :: B03_07, B03_11, B04_08, B04_12, B05_06, B05_07, B05_13, B06_08
    real :: B06_14, B07_08, B07_15, B08_16, B09_10, B09_11, B09_13, B10_12
    real :: B10_14, B11_12, B11_15, B12_16, B13_14, B13_15, B14_16, B15_16

    common U, w, w0, phi, w1, w2, w3, w4, J1, J2, J3, E01, E02, E03, E04, &
    E05, E06, E07, E08, E09, E10, E11, E12, E13, E14, E15, E16

    theta = w*t + phi
    B01_02 = (E01 - E02)*t + theta
    B01_03 = (E01 - E03)*t + theta
    B01_05 = (E01 - E05)*t + theta
    B01_09 = (E01 - E09)*t + theta

    dX01 = (U/2)*( (-X02*sin(B01_02))  + (-Y02*cos(B01_02))  &
			    +  (-X03*sin(B01_03))  + (-Y03*cos(B01_03))  &
			    +  (-X05*sin(B01_05))  + (-Y05*cos(B01_05))  &
			    +  (-X09*sin(B01_09))  + (-Y09*cos(B01_09))  )
end function

real function dX02(t, X01, Y01, X04, Y04, X06, Y06, X10, Y10)

    real :: U, theta, w, w0, phi,  w1, w2, w3, w4, J1, J2, J3
    real :: E01, E02, E03, E04, E05, E06, E07, E08, E09, E10, E11, E12, E13, E14, E15, E16
    real :: t, X01, Y01, X04, Y04, X06, Y06, X10, Y10

    real :: B01_02, B01_03, B01_05, B01_09, B02_04, B02_06, B02_10, B03_04
    real :: B03_07, B03_11, B04_08, B04_12, B05_06, B05_07, B05_13, B06_08
    real :: B06_14, B07_08, B07_15, B08_16, B09_10, B09_11, B09_13, B10_12
    real :: B10_14, B11_12, B11_15, B12_16, B13_14, B13_15, B14_16, B15_16

    common U, w, w0, phi, w1, w2, w3, w4, J1, J2, J3, E01, E02, E03, E04, &
    E05, E06, E07, E08, E09, E10, E11, E12, E13, E14, E15, E16

    theta = w*t + phi
    B01_02 = (E01 - E02)*t + theta
    B02_04 = (E02 - E04)*t + theta
    B02_06 = (E02 - E06)*t + theta
    B02_10 = (E02 - E10)*t + theta

    dX02 = (U/2)*( ( X01*sin(B01_02))  + (-Y01*cos(B01_02))  &
			    +  (-X04*sin(B02_04))  + (-Y04*cos(B02_04))  &
			    +  (-X06*sin(B02_06))  + (-Y06*cos(B02_06))  &
			    +  (-X10*sin(B02_10))  + (-Y10*cos(B02_10))  )
end function

real function dX03(t, X01, Y01, X04, Y04, X07, Y07, X11, Y11)

    real :: U, theta, w, w0, phi,  w1, w2, w3, w4, J1, J2, J3
    real :: E01, E02, E03, E04, E05, E06, E07, E08, E09, E10, E11, E12, E13, E14, E15, E16
    real :: t, X01, Y01, X04, Y04, X07, Y07, X11, Y11

    real :: B01_02, B01_03, B01_05, B01_09, B02_04, B02_06, B02_10, B03_04
    real :: B03_07, B03_11, B04_08, B04_12, B05_06, B05_07, B05_13, B06_08
    real :: B06_14, B07_08, B07_15, B08_16, B09_10, B09_11, B09_13, B10_12
    real :: B10_14, B11_12, B11_15, B12_16, B13_14, B13_15, B14_16, B15_16

    common U, w, w0, phi, w1, w2, w3, w4, J1, J2, J3, E01, E02, E03, E04, &
    E05, E06, E07, E08, E09, E10, E11, E12, E13, E14, E15, E16

    theta = w*t + phi
    B01_03 = (E01 - E03)*t + theta
    B03_04 = (E03 - E04)*t + theta
    B03_07 = (E03 - E07)*t + theta
    B03_11 = (E03 - E11)*t + theta

    dX03 = (U/2)*( ( X01*sin(B01_03))  + (-Y01*cos(B01_03))  &
			    +  (-X04*sin(B03_04))  + (-Y04*cos(B03_04))  &
			    +  (-X07*sin(B03_07))  + (-Y07*cos(B03_07))  &
			    +  (-X11*sin(B03_11))  + (-Y11*cos(B03_11)) )
end function

real function dX04(t, X02, Y02, X03, Y03, X08, Y08, X12, Y12)

    real :: U, theta, w, w0, phi,  w1, w2, w3, w4, J1, J2, J3
    real :: E01, E02, E03, E04, E05, E06, E07, E08, E09, E10, E11, E12, E13, E14, E15, E16
    real :: t, X02, Y02, X03, Y03, X08, Y08, X12, Y12

    real :: B01_02, B01_03, B01_05, B01_09, B02_04, B02_06, B02_10, B03_04
    real :: B03_07, B03_11, B04_08, B04_12, B05_06, B05_07, B05_13, B06_08
    real :: B06_14, B07_08, B07_15, B08_16, B09_10, B09_11, B09_13, B10_12
    real :: B10_14, B11_12, B11_15, B12_16, B13_14, B13_15, B14_16, B15_16

    common U, w, w0, phi, w1, w2, w3, w4, J1, J2, J3, E01, E02, E03, E04, &
    E05, E06, E07, E08, E09, E10, E11, E12, E13, E14, E15, E16

    theta = w*t + phi
    B02_04 = (E02 - E04)*t + theta
    B03_04 = (E03 - E04)*t + theta
    B04_08 = (E04 - E08)*t + theta
    B04_12 = (E04 - E12)*t + theta

    dX04 = (U/2)*( ( X02*sin(B02_04))  + (-Y02*cos(B02_04))  &
			    +  ( X03*sin(B03_04))  + (-Y03*cos(B03_04))  &
			    +  (-X08*sin(B04_08))  + (-Y08*cos(B04_08))  &
			    +  (-X12*sin(B04_12))  + (-Y12*cos(B04_12)) )
end function

real function dX05(t, X01, Y01, X06, Y06, X07, Y07, X13, Y13)

    real :: U, theta, w, w0, phi,  w1, w2, w3, w4, J1, J2, J3
    real :: E01, E02, E03, E04, E05, E06, E07, E08, E09, E10, E11, E12, E13, E14, E15, E16
    real :: t, X01, Y01, X06, Y06, X07, Y07, X13, Y13

    real :: B01_02, B01_03, B01_05, B01_09, B02_04, B02_06, B02_10, B03_04
    real :: B03_07, B03_11, B04_08, B04_12, B05_06, B05_07, B05_13, B06_08
    real :: B06_14, B07_08, B07_15, B08_16, B09_10, B09_11, B09_13, B10_12
    real :: B10_14, B11_12, B11_15, B12_16, B13_14, B13_15, B14_16, B15_16

    common U, w, w0, phi, w1, w2, w3, w4, J1, J2, J3, E01, E02, E03, E04, &
    E05, E06, E07, E08, E09, E10, E11, E12, E13, E14, E15, E16



    theta = w*t + phi
    B01_05 = (E01 - E05)*t + theta
    B05_06 = (E05 - E06)*t + theta
    B05_07 = (E05 - E07)*t + theta
    B05_13 = (E05 - E13)*t + theta


    dX05 = (U/2)*( ( X01*sin(B01_05))  + (-Y01*cos(B01_05))  &
			    +  (-X06*sin(B05_06))  + (-Y06*cos(B05_06))  &
			    +  (-X07*sin(B05_07))  + (-Y07*cos(B05_07))  &
			    +  (-X13*sin(B05_13))  + (-Y13*cos(B05_13)) )
end function

real function dX06(t, X02, Y02, X05, Y05, X08, Y08, X14, Y14)

    real :: U, theta, w, w0, phi,  w1, w2, w3, w4, J1, J2, J3
    real :: E01, E02, E03, E04, E05, E06, E07, E08, E09, E10, E11, E12, E13, E14, E15, E16
    real :: t, X02, Y02, X05, Y05, X08, Y08, X14, Y14

    real :: B01_02, B01_03, B01_05, B01_09, B02_04, B02_06, B02_10, B03_04
    real :: B03_07, B03_11, B04_08, B04_12, B05_06, B05_07, B05_13, B06_08
    real :: B06_14, B07_08, B07_15, B08_16, B09_10, B09_11, B09_13, B10_12
    real :: B10_14, B11_12, B11_15, B12_16, B13_14, B13_15, B14_16, B15_16

    common U, w, w0, phi, w1, w2, w3, w4, J1, J2, J3, E01, E02, E03, E04, &
    E05, E06, E07, E08, E09, E10, E11, E12, E13, E14, E15, E16

    theta = w*t + phi
    B02_06 = (E02 - E06)*t + theta
    B05_06 = (E05 - E06)*t + theta
    B06_08 = (E06 - E08)*t + theta
    B06_14 = (E06 - E14)*t + theta

    dX06 = (U/2)*( ( X02*sin(B02_06))  + (-Y02*cos(B02_06))  &
			+  ( X05*sin(B05_06))  + (-Y05*cos(B05_06))  &
			+  (-X08*sin(B06_08))  + (-Y08*cos(B06_08))  &
			+  (-X14*sin(B06_14)) + (-Y14*cos(B06_14)) )
end function

real function dX07(t, X03, Y03, X05, Y05, X08, Y08, X15, Y15)

    real :: U, theta, w, w0, phi,  w1, w2, w3, w4, J1, J2, J3
    real :: E01, E02, E03, E04, E05, E06, E07, E08, E09, E10, E11, E12, E13, E14, E15, E16
    real :: t, X03, Y03, X05, Y05, X08, Y08, X15, Y15

    real :: B01_02, B01_03, B01_05, B01_09, B02_04, B02_06, B02_10, B03_04
    real :: B03_07, B03_11, B04_08, B04_12, B05_06, B05_07, B05_13, B06_08
    real :: B06_14, B07_08, B07_15, B08_16, B09_10, B09_11, B09_13, B10_12
    real :: B10_14, B11_12, B11_15, B12_16, B13_14, B13_15, B14_16, B15_16

    common U, w, w0, phi, w1, w2, w3, w4, J1, J2, J3, E01, E02, E03, E04, &
    E05, E06, E07, E08, E09, E10, E11, E12, E13, E14, E15, E16

    theta = w*t + phi
    B03_07 = (E03 - E07)*t + theta
    B05_07 = (E05 - E07)*t + theta
    B07_08 = (E07 - E08)*t + theta
    B07_15 = (E07 - E15)*t + theta

    dX07 = (U/2)*( ( X03*sin(B03_07))  + (-Y03*cos(B03_07))  &
			+  ( X05*sin(B05_07))  + (-Y05*cos(B05_07))  &
			+  (-X08*sin(B07_08))  + (-Y08*cos(B07_08))  &
			+  (-X15*sin(B07_15)) + (-Y15*cos(B07_15)) )
end function

real function dX08(t, X04, Y04, X06, Y06, X07, Y07, X16, Y16)

    real :: U, theta, w, w0, phi,  w1, w2, w3, w4, J1, J2, J3
    real :: E01, E02, E03, E04, E05, E06, E07, E08, E09, E10, E11, E12, E13, E14, E15, E16
    real :: t, X04, Y04, X06, Y06, X07, Y07, X16, Y16

    real :: B01_02, B01_03, B01_05, B01_09, B02_04, B02_06, B02_10, B03_04
    real :: B03_07, B03_11, B04_08, B04_12, B05_06, B05_07, B05_13, B06_08
    real :: B06_14, B07_08, B07_15, B08_16, B09_10, B09_11, B09_13, B10_12
    real :: B10_14, B11_12, B11_15, B12_16, B13_14, B13_15, B14_16, B15_16

    common U, w, w0, phi, w1, w2, w3, w4, J1, J2, J3, E01, E02, E03, E04, &
    E05, E06, E07, E08, E09, E10, E11, E12, E13, E14, E15, E16

    theta = w*t + phi
    B04_08 = (E04 - E08)*t + theta
    B06_08 = (E06 - E08)*t + theta
    B07_08 = (E07 - E08)*t + theta
    B08_16 = (E08 - E16)*t + theta

    dX08 = (U/2)*( ( X04*sin(B04_08))  + (-Y04*cos(B04_08))  &
			+  ( X06*sin(B06_08))  + (-Y06*cos(B06_08))  &
			+  ( X07*sin(B07_08))  + (-Y07*cos(B07_08))  &
			+  (-X16*sin(B08_16)) + (-Y16*cos(B08_16)) )
end function

real function dX09(t, X01, Y01, X10, Y10, X11, Y11, X13, Y13)

    real :: U, theta, w, w0, phi,  w1, w2, w3, w4, J1, J2, J3
    real :: E01, E02, E03, E04, E05, E06, E07, E08, E09, E10, E11, E12, E13, E14, E15, E16
    real :: t, X01, Y01, X10, Y10, X11, Y11, X13, Y13

    real :: B01_02, B01_03, B01_05, B01_09, B02_04, B02_06, B02_10, B03_04
    real :: B03_07, B03_11, B04_08, B04_12, B05_06, B05_07, B05_13, B06_08
    real :: B06_14, B07_08, B07_15, B08_16, B09_10, B09_11, B09_13, B10_12
    real :: B10_14, B11_12, B11_15, B12_16, B13_14, B13_15, B14_16, B15_16

    common U, w, w0, phi, w1, w2, w3, w4, J1, J2, J3, E01, E02, E03, E04, &
    E05, E06, E07, E08, E09, E10, E11, E12, E13, E14, E15, E16

    theta = w*t + phi
    B01_09 = (E01 - E09)*t + theta
    B09_10 = (E09 - E10)*t + theta
    B09_11 = (E09 - E11)*t + theta
    B09_13 = (E09 - E13)*t + theta


    dX09 = (U/2)*( ( X01*sin(B01_09))  + (-Y01*cos(B01_09))  &
			+  (-X10*sin(B09_10)) + (-Y10*cos(B09_10)) &
			+  (-X11*sin(B09_11)) + (-Y11*cos(B09_11)) &
			+  (-X13*sin(B09_13)) + (-Y13*cos(B09_13)) )
end function

real function dX10(t, X02, Y02, X09, Y09, X12, Y12, X14, Y14)

    real :: U, theta, w, w0, phi,  w1, w2, w3, w4, J1, J2, J3
    real :: E01, E02, E03, E04, E05, E06, E07, E08, E09, E10, E11, E12, E13, E14, E15, E16
    real :: t, X02, Y02, X09, Y09, X12, Y12, X14, Y14

    real :: B01_02, B01_03, B01_05, B01_09, B02_04, B02_06, B02_10, B03_04
    real :: B03_07, B03_11, B04_08, B04_12, B05_06, B05_07, B05_13, B06_08
    real :: B06_14, B07_08, B07_15, B08_16, B09_10, B09_11, B09_13, B10_12
    real :: B10_14, B11_12, B11_15, B12_16, B13_14, B13_15, B14_16, B15_16

    common U, w, w0, phi, w1, w2, w3, w4, J1, J2, J3, E01, E02, E03, E04, &
    E05, E06, E07, E08, E09, E10, E11, E12, E13, E14, E15, E16

    theta = w*t + phi
    B02_10 = (E02 - E10)*t + theta
    B09_10 = (E09 - E10)*t + theta
    B10_12 = (E10 - E12)*t + theta
    B10_14 = (E10 - E14)*t + theta

    dX10 = (U/2)*( ( X02*sin(B02_10)) + (-Y02*cos(B02_10)) &
			+  ( X09*sin(B09_10)) + (-Y09*cos(B09_10)) &
			+  (-X12*sin(B10_12))+ (-Y12*cos(B10_12))&
			+  (-X14*sin(B10_14))+ (-Y14*cos(B10_14)))
end function

real function dX11(t, X03, Y03, X09, Y09, X12, Y12, X15, Y15)

    real :: U, theta, w, w0, phi,  w1, w2, w3, w4, J1, J2, J3
    real :: E01, E02, E03, E04, E05, E06, E07, E08, E09, E10, E11, E12, E13, E14, E15, E16
    real :: t, X03, Y03, X09, Y09, X12, Y12, X15, Y15

    real :: B01_02, B01_03, B01_05, B01_09, B02_04, B02_06, B02_10, B03_04
    real :: B03_07, B03_11, B04_08, B04_12, B05_06, B05_07, B05_13, B06_08
    real :: B06_14, B07_08, B07_15, B08_16, B09_10, B09_11, B09_13, B10_12
    real :: B10_14, B11_12, B11_15, B12_16, B13_14, B13_15, B14_16, B15_16

    common U, w, w0, phi, w1, w2, w3, w4, J1, J2, J3, E01, E02, E03, E04, &
    E05, E06, E07, E08, E09, E10, E11, E12, E13, E14, E15, E16

    theta = w*t + phi
    B03_11 = (E03 - E11)*t + theta
    B09_11 = (E09 - E11)*t + theta
    B11_12 = (E11 - E12)*t + theta
    B11_15 = (E11 - E15)*t + theta

    dX11 =  (U/2)*( ( X03*sin(B03_11)) + (-Y03*cos(B03_11)) &
			 +  ( X09*sin(B09_11)) + (-Y09*cos(B09_11)) &
			 +  (-X12*sin(B11_12))+ (-Y12*cos(B11_12))&
			 +  (-X15*sin(B11_15))+ (-Y15*cos(B11_15)))
end function

real function dX12(t, X04, Y04, X10, Y10, X11, Y11, X16, Y16)

    real :: U, theta, w, w0, phi,  w1, w2, w3, w4, J1, J2, J3
    real :: E01, E02, E03, E04, E05, E06, E07, E08, E09, E10, E11, E12, E13, E14, E15, E16
    real :: t, X04, Y04, X10, Y10, X11, Y11, X16, Y16

    real :: B01_02, B01_03, B01_05, B01_09, B02_04, B02_06, B02_10, B03_04
    real :: B03_07, B03_11, B04_08, B04_12, B05_06, B05_07, B05_13, B06_08
    real :: B06_14, B07_08, B07_15, B08_16, B09_10, B09_11, B09_13, B10_12
    real :: B10_14, B11_12, B11_15, B12_16, B13_14, B13_15, B14_16, B15_16

    common U, w, w0, phi, w1, w2, w3, w4, J1, J2, J3, E01, E02, E03, E04, &
    E05, E06, E07, E08, E09, E10, E11, E12, E13, E14, E15, E16

    theta = w*t + phi
    B04_12 = (E04 - E12)*t + theta
    B10_12 = (E10 - E12)*t + theta
    B11_12 = (E11 - E12)*t + theta
    B12_16 = (E12 - E16)*t + theta

    dX12  = (U/2)*( ( X04*sin(B04_12)) + (-Y04*cos(B04_12)) &
			 +  ( X10*sin(B10_12))+ (-Y10*cos(B10_12))&
			 +  ( X11*sin(B11_12))+ (-Y11*cos(B11_12))&
			 +  (-X16*sin(B12_16))+ (-Y16*cos(B12_16)))
end function

real function dX13(t, X05, Y05, X09, Y09, X14, Y14, X15, Y15)

    real :: U, theta, w, w0, phi,  w1, w2, w3, w4, J1, J2, J3
    real :: E01, E02, E03, E04, E05, E06, E07, E08, E09, E10, E11, E12, E13, E14, E15, E16
    real :: t, X05, Y05, X09, Y09, X14, Y14, X15, Y15

    real :: B01_02, B01_03, B01_05, B01_09, B02_04, B02_06, B02_10, B03_04
    real :: B03_07, B03_11, B04_08, B04_12, B05_06, B05_07, B05_13, B06_08
    real :: B06_14, B07_08, B07_15, B08_16, B09_10, B09_11, B09_13, B10_12
    real :: B10_14, B11_12, B11_15, B12_16, B13_14, B13_15, B14_16, B15_16

    common U, w, w0, phi, w1, w2, w3, w4, J1, J2, J3, E01, E02, E03, E04, &
    E05, E06, E07, E08, E09, E10, E11, E12, E13, E14, E15, E16

    theta = w*t + phi
    B05_13 = (E05 - E13)*t + theta
    B09_13 = (E09 - E13)*t + theta
    B13_14 = (E13 - E14)*t + theta
    B13_15 = (E13 - E15)*t + theta

    dX13 =  (U/2)*( ( X05*sin(B05_13)) + (-Y05*cos(B05_13)) &
			 +  ( X09*sin(B09_13)) + (-Y09*cos(B09_13)) &
			 +  (-X14*sin(B13_14))+ (-Y14*cos(B13_14))&
			 +  (-X15*sin(B13_15))+ (-Y15*cos(B13_15)))
end function

real function dX14(t, X06, Y06, X10, Y10, X13, Y13, X16, Y16)
    real :: U, theta, w, w0, phi,  w1, w2, w3, w4, J1, J2, J3
    real :: E01, E02, E03, E04, E05, E06, E07, E08, E09, E10, E11, E12, E13, E14, E15, E16
    real :: t, X06, Y06, X10, Y10, X13, Y13, X16, Y16

    real :: B01_02, B01_03, B01_05, B01_09, B02_04, B02_06, B02_10, B03_04
    real :: B03_07, B03_11, B04_08, B04_12, B05_06, B05_07, B05_13, B06_08
    real :: B06_14, B07_08, B07_15, B08_16, B09_10, B09_11, B09_13, B10_12
    real :: B10_14, B11_12, B11_15, B12_16, B13_14, B13_15, B14_16, B15_16

    common U, w, w0, phi, w1, w2, w3, w4, J1, J2, J3, E01, E02, E03, E04, &
    E05, E06, E07, E08, E09, E10, E11, E12, E13, E14, E15, E16

    theta = w*t + phi
    B06_14 = (E06 - E14)*t + theta
    B10_14 = (E10 - E14)*t + theta
    B13_14 = (E13 - E14)*t + theta
    B14_16 = (E14 - E16)*t + theta

    dX14 = (U/2)*( ( X06*sin(B06_14)) + (-Y06*cos(B06_14)) &
		    +  ( X10*sin(B10_14))+ (-Y10*cos(B10_14))&
		    +  ( X13*sin(B13_14))+ (-Y13*cos(B13_14))&
		    +  (-X16*sin(B14_16))+ (-Y16*cos(B14_16)))
end function

real function dX15(t, X07, Y07, X11, Y11, X13, Y13, X16, Y16)
    real :: U, theta, w, w0, phi,  w1, w2, w3, w4, J1, J2, J3
    real :: E01, E02, E03, E04, E05, E06, E07, E08, E09, E10, E11, E12, E13, E14, E15, E16
    real :: t, X07, Y07, X11, Y11, X13, Y13, X16, Y16

    real :: B01_02, B01_03, B01_05, B01_09, B02_04, B02_06, B02_10, B03_04
    real :: B03_07, B03_11, B04_08, B04_12, B05_06, B05_07, B05_13, B06_08
    real :: B06_14, B07_08, B07_15, B08_16, B09_10, B09_11, B09_13, B10_12
    real :: B10_14, B11_12, B11_15, B12_16, B13_14, B13_15, B14_16, B15_16

    common U, w, w0, phi, w1, w2, w3, w4, J1, J2, J3, E01, E02, E03, E04, &
    E05, E06, E07, E08, E09, E10, E11, E12, E13, E14, E15, E16

    theta = w*t + phi
    B07_15 = (E07 - E15)*t + theta
    B11_15 = (E11 - E15)*t + theta
    B13_15 = (E13 - E15)*t + theta
    B15_16 = (E15 - E16)*t + theta

    dX15 = (U/2)*( ( X07*sin(B07_15)) + (-Y07*cos(B07_15)) &
		    +  ( X11*sin(B11_15))+ (-Y11*cos(B11_15))&
		    +  ( X13*sin(B13_15))+ (-Y13*cos(B13_15))&
		    +  (-X16*sin(B15_16))+ (-Y16*cos(B15_16)))
end function

real function dX16(t, X08, Y08, X12, Y12, X14, Y14, X15, Y15)

    real :: U, theta, w, w0, phi,  w1, w2, w3, w4, J1, J2, J3
    real :: E01, E02, E03, E04, E05, E06, E07, E08, E09, E10, E11, E12, E13, E14, E15, E16
    real :: t, X08, Y08, X12, Y12, X14, Y14, X15, Y15

    real :: B01_02, B01_03, B01_05, B01_09, B02_04, B02_06, B02_10, B03_04
    real :: B03_07, B03_11, B04_08, B04_12, B05_06, B05_07, B05_13, B06_08
    real :: B06_14, B07_08, B07_15, B08_16, B09_10, B09_11, B09_13, B10_12
    real :: B10_14, B11_12, B11_15, B12_16, B13_14, B13_15, B14_16, B15_16

    common U, w, w0, phi, w1, w2, w3, w4, J1, J2, J3, E01, E02, E03, E04, &
    E05, E06, E07, E08, E09, E10, E11, E12, E13, E14, E15, E16

    theta = w*t + phi
    B08_16 = (E08 - E16)*t + theta
    B12_16 = (E12 - E16)*t + theta
    B14_16 = (E14 - E16)*t + theta
    B15_16 = (E15 - E16)*t + theta

    dX16 = (U/2)*( ( X08*sin(B08_16)) + (-Y08*cos(B08_16)) &
		    +  ( X12*sin(B12_16))+ (-Y12*cos(B12_16))&
		    +  ( X14*sin(B14_16))+ (-Y14*cos(B14_16))&
		    +  ( X15*sin(B15_16))+ (-Y15*cos(B15_16)))
end function

! For DY

real function dY01(t, X02, Y02, X03, Y03, X05, Y05, X09, Y09)

    real :: U, theta, w, w0, phi,  w1, w2, w3, w4, J1, J2, J3
    real :: E01, E02, E03, E04, E05, E06, E07, E08, E09, E10, E11, E12, E13, E14, E15, E16
    real :: t, X02, Y02, X03, Y03, X05, Y05, X09, Y09

    real :: B01_02, B01_03, B01_05, B01_09, B02_04, B02_06, B02_10, B03_04
    real :: B03_07, B03_11, B04_08, B04_12, B05_06, B05_07, B05_13, B06_08
    real :: B06_14, B07_08, B07_15, B08_16, B09_10, B09_11, B09_13, B10_12
    real :: B10_14, B11_12, B11_15, B12_16, B13_14, B13_15, B14_16, B15_16

    common U, w, w0, phi, w1, w2, w3, w4, J1, J2, J3, E01, E02, E03, E04, &
    E05, E06, E07, E08, E09, E10, E11, E12, E13, E14, E15, E16

    theta = w*t + phi
    B01_02 = (E01 - E02)*t + theta
    B01_03 = (E01 - E03)*t + theta
    B01_05 = (E01 - E05)*t + theta
    B01_09 = (E01 - E09)*t + theta

    dY01 = (U/2)*( ( X02*cos(B01_02))  + (-Y02*sin(B01_02))  &
			+  ( X03*cos(B01_03))  + (-Y03*sin(B01_03))  &
			+  ( X05*cos(B01_05))  + (-Y05*sin(B01_05))  &
			+  ( X09*cos(B01_09))  + (-Y09*sin(B01_09))  )

end function

real function dY02(t, X01, Y01, X04, Y04, X06, Y06, X10, Y10)

    real :: U, theta, w, w0, phi,  w1, w2, w3, w4, J1, J2, J3
    real :: E01, E02, E03, E04, E05, E06, E07, E08, E09, E10, E11, E12, E13, E14, E15, E16
    real :: t, X01, Y01, X04, Y04, X06, Y06, X10, Y10

    real :: B01_02, B01_03, B01_05, B01_09, B02_04, B02_06, B02_10, B03_04
    real :: B03_07, B03_11, B04_08, B04_12, B05_06, B05_07, B05_13, B06_08
    real :: B06_14, B07_08, B07_15, B08_16, B09_10, B09_11, B09_13, B10_12
    real :: B10_14, B11_12, B11_15, B12_16, B13_14, B13_15, B14_16, B15_16

    common U, w, w0, phi, w1, w2, w3, w4, J1, J2, J3, E01, E02, E03, E04, &
    E05, E06, E07, E08, E09, E10, E11, E12, E13, E14, E15, E16

    theta = w*t + phi
    B01_02 = (E01 - E02)*t + theta
    B02_04 = (E02 - E04)*t + theta
    B02_06 = (E02 - E06)*t + theta
    B02_10 = (E02 - E10)*t + theta

    dY02 = (U/2)*( ( X01*cos(B01_02))  + ( Y01*sin(B01_02))  &
			+  ( X04*cos(B02_04))  + (-Y04*sin(B02_04))  &
			+  ( X06*cos(B02_06))  + (-Y06*sin(B02_06))  &
			+  ( X10*cos(B02_10)) + (-Y10*sin(B02_10)) )



end function

real function dY03(t, X01, Y01, X04, Y04, X07, Y07, X11, Y11)

    real :: U, theta, w, w0, phi,  w1, w2, w3, w4, J1, J2, J3
    real :: E01, E02, E03, E04, E05, E06, E07, E08, E09, E10, E11, E12, E13, E14, E15, E16
    real :: t, X01, Y01, X04, Y04, X07, Y07, X11, Y11

    real :: B01_02, B01_03, B01_05, B01_09, B02_04, B02_06, B02_10, B03_04
    real :: B03_07, B03_11, B04_08, B04_12, B05_06, B05_07, B05_13, B06_08
    real :: B06_14, B07_08, B07_15, B08_16, B09_10, B09_11, B09_13, B10_12
    real :: B10_14, B11_12, B11_15, B12_16, B13_14, B13_15, B14_16, B15_16

    common U, w, w0, phi, w1, w2, w3, w4, J1, J2, J3, E01, E02, E03, E04, &
    E05, E06, E07, E08, E09, E10, E11, E12, E13, E14, E15, E16

    theta = w*t + phi
    B01_03 = (E01 - E03)*t + theta
    B03_04 = (E03 - E04)*t + theta
    B03_07 = (E03 - E07)*t + theta
    B03_11 = (E03 - E11)*t + theta

    dY03 = (U/2)*( ( X01*cos(B01_03))  + ( Y01*sin(B01_03))  &
			+  ( X04*cos(B03_04))  + (-Y04*sin(B03_04))  &
			+  ( X07*cos(B03_07))  + (-Y07*sin(B03_07))  &
			+  ( X11*cos(B03_11)) + (-Y11*sin(B03_11)) )




end function

real function dY04(t, X02, Y02, X03, Y03, X08, Y08, X12, Y12)

    real :: U, theta, w, w0, phi,  w1, w2, w3, w4, J1, J2, J3
    real :: E01, E02, E03, E04, E05, E06, E07, E08, E09, E10, E11, E12, E13, E14, E15, E16
    real :: t, X02, Y02, X03, Y03, X08, Y08, X12, Y12

    real :: B01_02, B01_03, B01_05, B01_09, B02_04, B02_06, B02_10, B03_04
    real :: B03_07, B03_11, B04_08, B04_12, B05_06, B05_07, B05_13, B06_08
    real :: B06_14, B07_08, B07_15, B08_16, B09_10, B09_11, B09_13, B10_12
    real :: B10_14, B11_12, B11_15, B12_16, B13_14, B13_15, B14_16, B15_16

    common U, w, w0, phi, w1, w2, w3, w4, J1, J2, J3, E01, E02, E03, E04, &
    E05, E06, E07, E08, E09, E10, E11, E12, E13, E14, E15, E16

    theta = w*t + phi
    B02_04 = (E02 - E04)*t + theta
    B03_04 = (E03 - E04)*t + theta
    B04_08 = (E04 - E08)*t + theta
    B04_12 = (E04 - E12)*t + theta

    dY04 = (U/2)*( ( X02*cos(B02_04))  + ( Y02*sin(B02_04))  &
			    +  ( X03*cos(B03_04))  + ( Y03*sin(B03_04))  &
			    +  ( X08*cos(B04_08))  + (-Y08*sin(B04_08))  &
			    +  ( X12*cos(B04_12))  + (-Y12*sin(B04_12)) )




end function

real function dY05(t, X01, Y01, X06, Y06, X07, Y07, X13, Y13)

    real :: U, theta, w, w0, phi,  w1, w2, w3, w4, J1, J2, J3
    real :: E01, E02, E03, E04, E05, E06, E07, E08, E09, E10, E11, E12, E13, E14, E15, E16
    real :: t, X01, Y01, X06, Y06, X07, Y07, X13, Y13

    real :: B01_02, B01_03, B01_05, B01_09, B02_04, B02_06, B02_10, B03_04
    real :: B03_07, B03_11, B04_08, B04_12, B05_06, B05_07, B05_13, B06_08
    real :: B06_14, B07_08, B07_15, B08_16, B09_10, B09_11, B09_13, B10_12
    real :: B10_14, B11_12, B11_15, B12_16, B13_14, B13_15, B14_16, B15_16

    common U, w, w0, phi, w1, w2, w3, w4, J1, J2, J3, E01, E02, E03, E04, &
    E05, E06, E07, E08, E09, E10, E11, E12, E13, E14, E15, E16

    theta = w*t + phi
    B01_05 = (E01 - E05)*t + theta
    B05_06 = (E05 - E06)*t + theta
    B05_07 = (E05 - E07)*t + theta
    B05_13 = (E05 - E13)*t + theta

    dY05 = (U/2)*( ( X01*cos(B01_05))  + ( Y01*sin(B01_05))  &
			+  ( X06*cos(B05_06))  + (-Y06*sin(B05_06))  &
			+  ( X07*cos(B05_07))  + (-Y07*sin(B05_07))  &
			+  ( X13*cos(B05_13)) + (-Y13*sin(B05_13)) )




end function

real function dY06(t, X02, Y02, X05, Y05, X08, Y08, X14, Y14)

    real :: U, theta, w, w0, phi,  w1, w2, w3, w4, J1, J2, J3
    real :: E01, E02, E03, E04, E05, E06, E07, E08, E09, E10, E11, E12, E13, E14, E15, E16
    real :: t, X02, Y02, X05, Y05, X08, Y08, X14, Y14

    real :: B01_02, B01_03, B01_05, B01_09, B02_04, B02_06, B02_10, B03_04
    real :: B03_07, B03_11, B04_08, B04_12, B05_06, B05_07, B05_13, B06_08
    real :: B06_14, B07_08, B07_15, B08_16, B09_10, B09_11, B09_13, B10_12
    real :: B10_14, B11_12, B11_15, B12_16, B13_14, B13_15, B14_16, B15_16

    common U, w, w0, phi, w1, w2, w3, w4, J1, J2, J3, E01, E02, E03, E04, &
    E05, E06, E07, E08, E09, E10, E11, E12, E13, E14, E15, E16

    theta = w*t + phi
    B02_06 = (E02 - E06)*t + theta
    B05_06 = (E05 - E06)*t + theta
    B06_08 = (E06 - E08)*t + theta
    B06_14 = (E06 - E14)*t + theta

    dY06 = (U/2)*( ( X02*cos(B02_06))  + ( Y02*sin(B02_06))  &
			+  ( X05*cos(B05_06))  + ( Y05*sin(B05_06))  &
			+  ( X08*cos(B06_08))  + (-Y08*sin(B06_08))  &
			+  ( X14*cos(B06_14)) + (-Y14*sin(B06_14)) )



end function

real function dY07(t, X03, Y03, X05, Y05, X08, Y08, X15, Y15)

    real :: U, theta, w, w0, phi,  w1, w2, w3, w4, J1, J2, J3
    real :: E01, E02, E03, E04, E05, E06, E07, E08, E09, E10, E11, E12, E13, E14, E15, E16
    real :: t, X03, Y03, X05, Y05, X08, Y08, X15, Y15

    real :: B01_02, B01_03, B01_05, B01_09, B02_04, B02_06, B02_10, B03_04
    real :: B03_07, B03_11, B04_08, B04_12, B05_06, B05_07, B05_13, B06_08
    real :: B06_14, B07_08, B07_15, B08_16, B09_10, B09_11, B09_13, B10_12
    real :: B10_14, B11_12, B11_15, B12_16, B13_14, B13_15, B14_16, B15_16

    common U, w, w0, phi, w1, w2, w3, w4, J1, J2, J3, E01, E02, E03, E04, &
    E05, E06, E07, E08, E09, E10, E11, E12, E13, E14, E15, E16

    theta = w*t + phi
    B03_07 = (E03 - E07)*t + theta
    B05_07 = (E05 - E07)*t + theta
    B07_08 = (E07 - E08)*t + theta
    B07_15 = (E07 - E15)*t + theta

    dY07 = (U/2)*( ( X03*cos(B03_07))  + ( Y03*sin(B03_07))  &
			+  ( X05*cos(B05_07))  + ( Y05*sin(B05_07))  &
			+  ( X08*cos(B07_08))  + (-Y08*sin(B07_08))  &
			+  ( X15*cos(B07_15)) + (-Y15*sin(B07_15)) )



end function

real function dY08(t, X04, Y04, X06, Y06, X07, Y07, X16, Y16)

    real :: U, theta, w, w0, phi,  w1, w2, w3, w4, J1, J2, J3
    real :: E01, E02, E03, E04, E05, E06, E07, E08, E09, E10, E11, E12, E13, E14, E15, E16
    real :: t, X04, Y04, X06, Y06, X07, Y07, X16, Y16

    real :: B01_02, B01_03, B01_05, B01_09, B02_04, B02_06, B02_10, B03_04
    real :: B03_07, B03_11, B04_08, B04_12, B05_06, B05_07, B05_13, B06_08
    real :: B06_14, B07_08, B07_15, B08_16, B09_10, B09_11, B09_13, B10_12
    real :: B10_14, B11_12, B11_15, B12_16, B13_14, B13_15, B14_16, B15_16

    common U, w, w0, phi, w1, w2, w3, w4, J1, J2, J3, E01, E02, E03, E04, &
    E05, E06, E07, E08, E09, E10, E11, E12, E13, E14, E15, E16

    theta = w*t + phi
    B04_08 = (E04 - E08)*t + theta
    B06_08 = (E06 - E08)*t + theta
    B07_08 = (E07 - E08)*t + theta
    B08_16 = (E08 - E16)*t + theta

    dY08 = (U/2)*( ( X04*cos(B04_08))  + ( Y04*sin(B04_08))  &
			+  ( X06*cos(B06_08))  + ( Y06*sin(B06_08))  &
			+  ( X07*cos(B07_08))  + ( Y07*sin(B07_08))  &
			+  ( X16*cos(B08_16)) + (-Y16*sin(B08_16)) )



end function

real function dY09(t, X01, Y01, X10, Y10, X11, Y11, X13, Y13)

    real :: U, theta, w, w0, phi,  w1, w2, w3, w4, J1, J2, J3
    real :: E01, E02, E03, E04, E05, E06, E07, E08, E09, E10, E11, E12, E13, E14, E15, E16
    real :: t, X01, Y01, X10, Y10, X11, Y11, X13, Y13

    real :: B01_02, B01_03, B01_05, B01_09, B02_04, B02_06, B02_10, B03_04
    real :: B03_07, B03_11, B04_08, B04_12, B05_06, B05_07, B05_13, B06_08
    real :: B06_14, B07_08, B07_15, B08_16, B09_10, B09_11, B09_13, B10_12
    real :: B10_14, B11_12, B11_15, B12_16, B13_14, B13_15, B14_16, B15_16

    common U, w, w0, phi, w1, w2, w3, w4, J1, J2, J3, E01, E02, E03, E04, &
    E05, E06, E07, E08, E09, E10, E11, E12, E13, E14, E15, E16

    theta = w*t + phi
    B01_09 = (E01 - E09)*t + theta
    B09_10 = (E09 - E10)*t + theta
    B09_11 = (E09 - E11)*t + theta
    B09_13 = (E09 - E13)*t + theta

    dY09 = (U/2)*( ( X01*cos(B01_09))  + ( Y01*sin(B01_09))  &
			+  ( X10*cos(B09_10)) + (-Y10*sin(B09_10)) &
			+  ( X11*cos(B09_11)) + (-Y11*sin(B09_11)) &
			+  ( X13*cos(B09_13)) + (-Y13*sin(B09_13)) )



end function

real function dY10(t, X02, Y02, X09, Y09, X12, Y12, X14, Y14)

    real :: U, theta, w, w0, phi,  w1, w2, w3, w4, J1, J2, J3
    real :: E01, E02, E03, E04, E05, E06, E07, E08, E09, E10, E11, E12, E13, E14, E15, E16
    real :: t, X02, Y02, X09, Y09, X12, Y12, X14, Y14

    real :: B01_02, B01_03, B01_05, B01_09, B02_04, B02_06, B02_10, B03_04
    real :: B03_07, B03_11, B04_08, B04_12, B05_06, B05_07, B05_13, B06_08
    real :: B06_14, B07_08, B07_15, B08_16, B09_10, B09_11, B09_13, B10_12
    real :: B10_14, B11_12, B11_15, B12_16, B13_14, B13_15, B14_16, B15_16

    common U, w, w0, phi, w1, w2, w3, w4, J1, J2, J3, E01, E02, E03, E04, &
    E05, E06, E07, E08, E09, E10, E11, E12, E13, E14, E15, E16

    theta = w*t + phi
    B02_10 = (E02 - E10)*t + theta
    B09_10 = (E09 - E10)*t + theta
    B10_12 = (E10 - E12)*t + theta
    B10_14 = (E10 - E14)*t + theta

    dY10 =  (U/2)*( ( X02*cos(B02_10)) + ( Y02*sin(B02_10)) &
			+  ( X09*cos(B09_10)) + ( Y09*sin(B09_10)) &
			+  ( X12*cos(B10_12))+ (-Y12*sin(B10_12))&
			+  ( X14*cos(B10_14))+ (-Y14*sin(B10_14)))
end function

real function dY11(t, X03, Y03, X09, Y09, X12, Y12, X15, Y15)

    real :: U, theta, w, w0, phi,  w1, w2, w3, w4, J1, J2, J3
    real :: E01, E02, E03, E04, E05, E06, E07, E08, E09, E10, E11, E12, E13, E14, E15, E16
    real :: t, X03, Y03, X09, Y09, X12, Y12, X15, Y15

    real :: B01_02, B01_03, B01_05, B01_09, B02_04, B02_06, B02_10, B03_04
    real :: B03_07, B03_11, B04_08, B04_12, B05_06, B05_07, B05_13, B06_08
    real :: B06_14, B07_08, B07_15, B08_16, B09_10, B09_11, B09_13, B10_12
    real :: B10_14, B11_12, B11_15, B12_16, B13_14, B13_15, B14_16, B15_16

    common U, w, w0, phi, w1, w2, w3, w4, J1, J2, J3, E01, E02, E03, E04, &
    E05, E06, E07, E08, E09, E10, E11, E12, E13, E14, E15, E16

    theta = w*t + phi
    B03_11 = (E03 - E11)*t + theta
    B09_11 = (E09 - E11)*t + theta
    B11_12 = (E11 - E12)*t + theta
    B11_15 = (E11 - E15)*t + theta

    dY11 =  (U/2)*( ( X03*cos(B03_11)) + ( Y03*sin(B03_11)) &
			 +  ( X09*cos(B09_11)) + ( Y09*sin(B09_11)) &
			 +  ( X12*cos(B11_12))+ (-Y12*sin(B11_12))&
			 +  ( X15*cos(B11_15))+ (-Y15*sin(B11_15)))



end function

real function dY12(t, X04, Y04, X10, Y10, X11, Y11, X16, Y16)

    real :: U, theta, w, w0, phi,  w1, w2, w3, w4, J1, J2, J3
    real :: E01, E02, E03, E04, E05, E06, E07, E08, E09, E10, E11, E12, E13, E14, E15, E16
    real :: t, X04, Y04, X10, Y10, X11, Y11, X16, Y16

    real :: B01_02, B01_03, B01_05, B01_09, B02_04, B02_06, B02_10, B03_04
    real :: B03_07, B03_11, B04_08, B04_12, B05_06, B05_07, B05_13, B06_08
    real :: B06_14, B07_08, B07_15, B08_16, B09_10, B09_11, B09_13, B10_12
    real :: B10_14, B11_12, B11_15, B12_16, B13_14, B13_15, B14_16, B15_16

    common U, w, w0, phi, w1, w2, w3, w4, J1, J2, J3, E01, E02, E03, E04, &
    E05, E06, E07, E08, E09, E10, E11, E12, E13, E14, E15, E16

    theta = w*t + phi
    B04_12 = (E04 - E12)*t + theta
    B10_12 = (E10 - E12)*t + theta
    B11_12 = (E11 - E12)*t + theta
    B12_16 = (E12 - E16)*t + theta

    dY12 =  (U/2)*( ( X04*cos(B04_12)) + ( Y04*sin(B04_12)) &
			 +  ( X10*cos(B10_12))+ ( Y10*sin(B10_12))&
			 +  ( X11*cos(B11_12))+ ( Y11*sin(B11_12))&
			 +  ( X16*cos(B12_16))+ (-Y16*sin(B12_16)))



end function

real function dY13(t, X05, Y05, X09, Y09, X14, Y14, X15, Y15)

    real :: U, theta, w, w0, phi,  w1, w2, w3, w4, J1, J2, J3
    real :: E01, E02, E03, E04, E05, E06, E07, E08, E09, E10, E11, E12, E13, E14, E15, E16
    real :: t, X05, Y05, X09, Y09, X14, Y14, X15, Y15

    real :: B01_02, B01_03, B01_05, B01_09, B02_04, B02_06, B02_10, B03_04
    real :: B03_07, B03_11, B04_08, B04_12, B05_06, B05_07, B05_13, B06_08
    real :: B06_14, B07_08, B07_15, B08_16, B09_10, B09_11, B09_13, B10_12
    real :: B10_14, B11_12, B11_15, B12_16, B13_14, B13_15, B14_16, B15_16

    common U, w, w0, phi, w1, w2, w3, w4, J1, J2, J3, E01, E02, E03, E04, &
    E05, E06, E07, E08, E09, E10, E11, E12, E13, E14, E15, E16

    theta = w*t + phi
    B05_13 = (E05 - E13)*t + theta
    B09_13 = (E09 - E13)*t + theta
    B13_14 = (E13 - E14)*t + theta
    B13_15 = (E13 - E15)*t + theta

    dY13 =  (U/2)*( ( X05*cos(B05_13)) + ( Y05*sin(B05_13)) &
			 +  ( X09*cos(B09_13)) + ( Y09*sin(B09_13)) &
			 +  ( X14*cos(B13_14))+ (-Y14*sin(B13_14))&
			 +  ( X15*cos(B13_15))+ (-Y15*sin(B13_15)))



end function

real function dY14(t, X06, Y06, X10, Y10, X13, Y13, X16, Y16)

    real :: U, theta, w, w0, phi,  w1, w2, w3, w4, J1, J2, J3
    real :: E01, E02, E03, E04, E05, E06, E07, E08, E09, E10, E11, E12, E13, E14, E15, E16
    real :: t, X06, Y06, X10, Y10, X13, Y13, X16, Y16

    real :: B01_02, B01_03, B01_05, B01_09, B02_04, B02_06, B02_10, B03_04
    real :: B03_07, B03_11, B04_08, B04_12, B05_06, B05_07, B05_13, B06_08
    real :: B06_14, B07_08, B07_15, B08_16, B09_10, B09_11, B09_13, B10_12
    real :: B10_14, B11_12, B11_15, B12_16, B13_14, B13_15, B14_16, B15_16

    common U, w, w0, phi, w1, w2, w3, w4, J1, J2, J3, E01, E02, E03, E04, &
    E05, E06, E07, E08, E09, E10, E11, E12, E13, E14, E15, E16

    theta = w*t + phi
    B06_14 = (E06 - E14)*t + theta
    B10_14 = (E10 - E14)*t + theta
    B13_14 = (E13 - E14)*t + theta
    B14_16 = (E14 - E16)*t + theta

    dY14 = (U/2)*( ( X06*cos(B06_14)) + ( Y06*sin(B06_14)) &
		    +  ( X10*cos(B10_14))+ ( Y10*sin(B10_14))&
		    +  ( X13*cos(B13_14))+ ( Y13*sin(B13_14))&
		    +  ( X16*cos(B14_16))+ (-Y16*sin(B14_16)))



end function

real function dY15(t, X07, Y07, X11, Y11, X13, Y13, X16, Y16)

    real :: U, theta, w, w0, phi,  w1, w2, w3, w4, J1, J2, J3
    real :: E01, E02, E03, E04, E05, E06, E07, E08, E09, E10, E11, E12, E13, E14, E15, E16
    real :: t, X07, Y07, X11, Y11, X13, Y13, X16, Y16

    real :: B01_02, B01_03, B01_05, B01_09, B02_04, B02_06, B02_10, B03_04
    real :: B03_07, B03_11, B04_08, B04_12, B05_06, B05_07, B05_13, B06_08
    real :: B06_14, B07_08, B07_15, B08_16, B09_10, B09_11, B09_13, B10_12
    real :: B10_14, B11_12, B11_15, B12_16, B13_14, B13_15, B14_16, B15_16

    common U, w, w0, phi, w1, w2, w3, w4, J1, J2, J3, E01, E02, E03, E04, &
    E05, E06, E07, E08, E09, E10, E11, E12, E13, E14, E15, E16

    theta = w*t + phi
    B07_15 = (E07 - E15)*t + theta
    B11_15 = (E11 - E15)*t + theta
    B13_15 = (E13 - E15)*t + theta
    B15_16 = (E15 - E16)*t + theta

    dY15 = (U/2)*( ( X07*cos(B07_15)) + ( Y07*sin(B07_15)) &
		    +  ( X11*cos(B11_15))+ ( Y11*sin(B11_15))&
		    +  ( X13*cos(B13_15))+ ( Y13*sin(B13_15))&
		    +  ( X16*cos(B15_16))+ (-Y16*sin(B15_16)))


end function

real function dY16(t, X08, Y08, X12, Y12, X14, Y14, X15, Y15)

    real :: U, theta, w, w0, phi,  w1, w2, w3, w4, J1, J2, J3
    real :: E01, E02, E03, E04, E05, E06, E07, E08, E09, E10, E11, E12, E13, E14, E15, E16
    real :: t, X08, Y08, X12, Y12, X14, Y14, X15, Y15

    real :: B01_02, B01_03, B01_05, B01_09, B02_04, B02_06, B02_10, B03_04
    real :: B03_07, B03_11, B04_08, B04_12, B05_06, B05_07, B05_13, B06_08
    real :: B06_14, B07_08, B07_15, B08_16, B09_10, B09_11, B09_13, B10_12
    real :: B10_14, B11_12, B11_15, B12_16, B13_14, B13_15, B14_16, B15_16

    common U, w, w0, phi, w1, w2, w3, w4, J1, J2, J3, E01, E02, E03, E04, &
    E05, E06, E07, E08, E09, E10, E11, E12, E13, E14, E15, E16

    theta = w*t + phi
    B08_16 = (E08 - E16)*t + theta
    B12_16 = (E12 - E16)*t + theta
    B14_16 = (E14 - E16)*t + theta
    B15_16 = (E15 - E16)*t + theta

    dY16 = (U/2)*( ( X08*cos(B08_16)) + ( Y08*sin(B08_16)) &
		    +  ( X12*cos(B12_16))+ ( Y12*sin(B12_16))&
		    +  ( X14*cos(B14_16))+ ( Y14*sin(B14_16))&
		    +  ( X15*cos(B15_16))+ ( Y15*sin(B15_16)))


end function

