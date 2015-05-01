program ciclo

use loading
implicit none
character(10)::Strings,str 

teste_var = 0
var_dois = 0

if (teste_var == 0 .or. var_dois == 0) then
        write(*,*) 'Passou no ciclo'
end if

! Teste com strings

Strings = 'HexTretas'
str = ' ' 
write(*,*) Strings
str = read(3,*) Strings
write(*,*) str
end program
