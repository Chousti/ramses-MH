MODULE imf_module
   ! Module for dealing with IMF calculations
   use amr_parameters, only: dp
   implicit none

   public :: sample_IMF_pop3
CONTAINS

FUNCTION sample_IMF_pop3() result(mass)
   !! Samples pop III IMF
   use pm_parameters, only: p3_mchar
   use pm_commons, ONLY: localseed
   use random
   implicit none
   real(dp):: mass
   real(dp):: total_mass
   real(dp):: m_min, m_max, d_log10m
   real(dp):: total_probability, local_probability
   real(dp)::RandNum
   integer:: n_samples = 1000
   integer:: i

   ! Min and Max of the IMF
   m_min = 1.d0
   m_max = 1000.d0

   ! Get the total_probability
   d_log10m = (LOG10(m_max) - LOG10(m_min)) / n_samples

   ! Get the total_probability
   total_probability = 0.d0
   do i=1,n_samples
      mass = 10.d0 ** (LOG10(m_min) + d_log10m * real(i-1,kind=dp))
      total_probability = total_probability + (mass**-2.3d0) * exp(-(p3_mchar/mass)**1.6d0)
   end do

   ! Loop again to the mass
   local_probability = 0.0
   mass = 0.0
   call random_number(RandNum)
   i = 1
   do while (local_probability.lt.RandNum)
      mass = 10.d0 ** (LOG10(m_min) + d_log10m * real(i-1,kind=dp))
      local_probability = local_probability + ((mass**-2.3d0) * exp(-(p3_mchar/mass)**1.6d0))/total_probability
      i = i + 1
   end do

END FUNCTION sample_IMF_pop3


END MODULE imf_module