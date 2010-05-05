#
# This file and all modifications and additions to the pristine
# package are under the same license as the package itself.
#

# norootforbuild

Name:           engrid
Version:	1.2rc1
Release:	1
Summary:	open-source mesh generator for CFD
Group:		Applications/Engineering
License:	GPL
Url:		http://engits.eu/engrid
Requires:	libqt4 vtk
BuildRequires:  libqt4-devel vtk netgen
Source:		http://files.engits.eu/engrid-1.2rc1.tar.gz
BuildRoot:      %{_tmppath}/%{name}-%{version}-build

%description

%prep
%setup

%build
cd src
export CXXFLAGS="$RPM_OPT_FLAGS"
export VTKINCDIR=%_includedir/vtk
export VTKLIBDIR=%_libdir
scripts/build-all.sh

%install
cp src/engrid /usr/bin

%clean
rm -rf $RPM_BUILD_ROOT

%post
%postun

%files
/usr/bin/engrid

%changelog
* Wed May 05 2010 ogloth@engits.com
- first attempt to create RPM package for enGrid


