class Rastair < Formula
  desc "A Rust-based command-line tool for genomic data processing"
  homepage "https://docs.rastair.com/"
  url "FILL_THIS_IN"
  sha256 "FILL_THIS_IN"
  license "AGPL"
  version "2.0.0-rc.1"

  # Declare build and library dependencies
  depends_on "rust" => :build
  depends_on "cmake" => :build
  depends_on "htslib"
  depends_on "libdeflate"

  def install
    # Set environment variables so your Rust tool finds the Homebrew libraries
    ENV.prepend_path "PKG_CONFIG_PATH", Formula["htslib"].opt_lib/"pkgconfig"
    ENV.prepend_path "PKG_CONFIG_PATH", Formula["libdeflate"].opt_lib/"pkgconfig"

    # Build and install the Rust executable
    system "cargo", "install", *std_cargo_args

    # Generate and install shell completions
    bash_completion_file = buildpath/"rastair.bash"
    File.write(bash_completion_file, `#{bin}/rastair internal shell-completions bash`)
    bash_completion.install bash_completion_file

    zsh_completion_file = buildpath/"rastair.zsh"
    File.write(zsh_completion_file, `#{bin}/rastair internal shell-completions zsh`)
    zsh_completion.install zsh_completion_file

    fish_completion_file = buildpath/"rastair.fish"
    File.write(fish_completion_file, `#{bin}/rastair internal shell-completions fish`)
    fish_completion.install fish_completion_file
  end

  test do
    # Example test to verify the binary works
    assert_match version.to_s, shell_output("#{bin}/rastair --version")
  end
end
