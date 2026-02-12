class Rastair < Formula
  desc "Rust-based command-line tool for genomic data processing"
  homepage "https://www.rastair.com/"
  license :cannot_represent

  if OS.mac?
    if Hardware::CPU.arm?
      url "https://s3.eu-west-2.amazonaws.com/com.rastair.releases/build/release-v2.0.0/rastair-v2.0.0-aarch64-apple-darwin.zip"
      sha256 "359c3c9432ec6e4d537aaf2bb86740483bf844b401884b132c66b63be965e826"
    else
      url "https://s3.eu-west-2.amazonaws.com/com.rastair.releases/build/release-v2.0.0/rastair-v2.0.0-x86_64-apple-darwin.zip"
      sha256 "814f64fd43da118eda71ece56ab39980ed3d5c88aa2b53a11b975e11d49a5d8f"
    end
  elsif OS.linux?
    url "https://s3.eu-west-2.amazonaws.com/com.rastair.releases/build/release-v2.0.0/rastair-v2.0.0-x86_64-unknown-linux-gnu.tar.gz"
    sha256 "8e74678a719e7c1ffb7ddc8a0fc9091878dd946da8610652b9151b7f25c9dcbe"
  end

  head do
    url "https://bitbucket.org/bsblabludwig/rastair.git", branch: "main"
    depends_on "cmake" => :build
    depends_on "htslib" => :build
    depends_on "libdeflate" => :build
    depends_on "rust" => :build
  end

  def install
    if build.head?
      # Set environment variables so Rust tool finds the Homebrew libraries
      ENV.prepend_path "PKG_CONFIG_PATH", Formula["htslib"].opt_lib/"pkgconfig"
      ENV.prepend_path "PKG_CONFIG_PATH", Formula["libdeflate"].opt_lib/"pkgconfig"

      system "cargo", "install", *std_cargo_args

      pkgshare.install "scripts/mbias.R", "scripts/QC_report.Rmd"
    else
      bin.install "rastair"
      pkgshare.install "mbias.R", "QC_report.Rmd"
    end

    # Generate and install shell completions
    bash_completion_file = buildpath/"rastair.bash"
    File.write(bash_completion_file, Utils.safe_popen_read(bin/"rastair", "internal", "shell-completions", "bash"))
    bash_completion.install bash_completion_file

    zsh_completion_file = buildpath/"rastair.zsh"
    File.write(zsh_completion_file, Utils.safe_popen_read(bin/"rastair", "internal", "shell-completions", "zsh"))
    zsh_completion.install zsh_completion_file

    fish_completion_file = buildpath/"rastair.fish"
    File.write(fish_completion_file, Utils.safe_popen_read(bin/"rastair", "internal", "shell-completions", "fish"))
    fish_completion.install fish_completion_file
  end

  test do
    assert_match version.to_s, shell_output("#{bin}/rastair --version")
  end
end
