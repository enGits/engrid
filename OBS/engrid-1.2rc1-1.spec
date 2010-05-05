#
# This file and all modifications and additions to the pristine
# package are under the same license as the package itself.
#

# norootforbuild

Name:           engrid-1.2rc1-1
Version:	1.2rc1
Release:	1
Summary:	open-source mesh generator for CFD
Group:		Applications/Engineering
License:	GPL
Url:		http://engits.eu/engrid
Requires:	libqt4 vtk-qt
BuildRequires:  libqt4-devel vtk-devel
Source:		http://files.engits.eu/engrid-1.2rc1.tar.gz
Patch:
BuildRoot:      %{_tmppath}/%{name}-%{version}-build

%description

%prep
%setup

%build

%install

%clean
rm -rf $RPM_BUILD_ROOT

%post
%postun

%files
%defattr(-,root,root)
%doc ChangeLog README COPYING

%changelog
* Wed May 05 2010 ogloth@ogloth-LAPTOP

