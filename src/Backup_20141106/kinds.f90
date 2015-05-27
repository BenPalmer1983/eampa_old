Module kinds
  Implicit None
  Integer, Parameter :: SingleReal = Selected_Real_Kind(6,37)        ! single real, 6 decimal precision, exponent range 37    
  Integer, Parameter :: DoubleReal = Selected_Real_Kind(15,307)      ! double real, 15 decimal precision, exponent range 307    
  Integer, Parameter :: QuadrupoleReal = Selected_Real_Kind(33,4931)  ! quadrupole real
  Integer, Parameter :: TinyInteger = Selected_Int_Kind(1)            ! tiny integer    1 byte
  Integer, Parameter :: SmallInteger = Selected_Int_Kind(4)           ! small integer    4 bytes -2E31 to 2E31-1
  Integer, Parameter :: StandardInteger = Selected_Int_Kind(8)        ! standard integer 8 bytes -2E63 to 2E63-1
  Integer, Parameter :: LongInteger = Selected_Int_Kind(12)           ! long integer
  Integer, Parameter :: VeryLongInteger = Selected_Int_Kind(32)       ! very long integer
  
  !example usage
  ! Integer(kind=StandardInteger) ::
  ! Real(kind=SingleReal) ::
  ! Real(kind=DoubleReal) ::
  
End Module kinds